# Figure 3: raw-source final audit (R, Windows)
# Target manuscript: Communications Biology | NOT a new GWAS or historical workflow reproduction
# Outputs are deposited in a new timestamped audit directory; manuscript/figures untouched.
# Inputs: two *original* GWAS .gz files and the archived MVP candidate plan + representative SNP table.
# This is a documented POST HOC reanalysis. Never tune inclusion criteria to force a historical match.

suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("CDD_HTN_ROOT", unset = ".")
root <- normalizePath(root, mustWork = TRUE)
htn_path <- file.path(root, '01_genetics/GWAS/HTN/GCST90044351_buildGRCh37.tsv.gz')
mvp_path <- file.path(root, '04_replication/HTN_BP_MVP/raw/MVP_R4.1000G_AGR.Systolic_Mean_INT.EUR.GIA.dbGaP.txt.gz')
repdir <- file.path(root, '04_replication/HTN_BP_MVP')
plan_path <- file.path(repdir, 'MVP_replication_plan_frozen_2026-09-13.tsv')
lead_path <- file.path(repdir, 'MVP_SBP_replication_harmonized_2026-09-13.tsv')
hist_path <- file.path(repdir, 'MVP_SBP_replication_FINAL_summary_2026-09-13.tsv')
inputs <- c(htn=htn_path, mvp=mvp_path, plan=plan_path, representatives=lead_path, historical=hist_path)
if (!all(file.exists(inputs))) stop('Missing input(s): ', paste(inputs[!file.exists(inputs)], collapse=' | '))

out <- file.path(root, '09_manuscript_audit', paste0('Figure3_R_raw_final_audit_', format(Sys.time(), '%Y%m%d_%H%M%S')))
if (dir.exists(out)) stop('Refusing to overwrite: ', out)
dir.create(out, recursive=TRUE, showWarnings=FALSE)
logfile <- file.path(out, '00_console_log.txt')
message('Output directory: ', out)
writeLines(c('R raw-GWAS audit. Historical Figure 3b procedures are not presumed recovered.',
             paste('Started', Sys.time()), paste('HTN:',htn_path),paste('MVP:',mvp_path)), logfile)
report <- function(...) { s <- paste(..., collapse=' '); cat(s, '\n'); cat(s, '\n', file=logfile, append=TRUE) }
save_tsv <- function(obj, name) fwrite(obj, file.path(out,name), sep='\t', na='NA')

# Locked audit choices (post hoc): GRCh37 +/- 1 Mb; unique rsID; same chromosome;
# direct A/C/G/T allele set only; valid beta, SE, P, EAF; |aligned EAF difference| <=0.10;
# main regional analysis removes palindromic A/T and C/G. Palindrome sensitivity is separate.
# Primary/representative SNPs are checked in a separate set, including frequency-resolvable palindromes.
plan <- fread(plan_path)
req_plan <- c('primary_locus_index','topSNP','topSNP_chr','topSNP_bp','candidate_region')
stopifnot(all(req_plan %in% names(plan)))
index <- plan[!is.na(primary_locus_index) & primary_locus_index == TRUE,
              .(region=candidate_region, chr=as.character(topSNP_chr), index_snp=topSNP, center=as.numeric(topSNP_bp))]
if (nrow(index)!=2L || !setequal(index$index_snp,c('rs1060240','rs1788799')) ||
    !setequal(index$chr,c('10','18'))) stop('Archived primary-index scheme differs; inspect plan, do not auto-pick variants.')
index[, `:=`(left=center-1000000, right=center+1000000)]
setorder(index,chr)
save_tsv(index,'01_locked_audit_windows.tsv')
report('AUDIT WINDOWS:', paste(index$region,index$chr,index$left,index$right,collapse=' ; '))

# The source gzip is sometimes described as concatenated gzip/BGZF. Read to actual EOF;
# if fewer than 1 million source lines have been scanned or an index SNP is absent, STOP.
# This prevents silently treating the first 65,280 uncompressed bytes as the whole file.
extract_stream <- function(path, mode, sought=NULL, chunk_size=50000L) {
  con <- gzfile(path,open='rt'); on.exit(close(con),add=TRUE)
  hdr <- readLines(con,n=1L,warn=FALSE)
  if (length(hdr)!=1L) stop(mode,': empty/unreadable gzip header')
  names_raw <- strsplit(hdr,'\t',fixed=TRUE)[[1]]
  expected <- if (mode=='HTN') c('chromosome','variant_id','base_pair_location','effect_allele','other_allele','beta','standard_error','p_value')
              else c('SNP_ID','chrom','ref','alt','ea','af','beta','sebeta','pval')
  if (!all(expected %in% names_raw)) stop(mode,': tab-delimited header not as expected; inspect original file manually.')
  report(mode, 'header verified; starting complete gzip scan ...')
  selected <- list(); batches<-0L; total<-0; inspected<-0; nmatch<-0L
  repeat {
    lines <- readLines(con,n=chunk_size,warn=FALSE)
    if (!length(lines)) break
    inspected <- inspected + length(lines)
    if (mode=='HTN') {
      chr10 <- startsWith(lines,'10\t'); chr18 <- startsWith(lines,'18\t')
      i <- which(chr10 | chr18)
      if (length(i)) {
        pos <- suppressWarnings(as.numeric(sub('^[^\t]*\t[^\t]*\t([^\t]+).*$', '\\1',lines[i],perl=TRUE)))
        loc <- ifelse(chr10[i], '10', '18'); q<-match(loc,index$chr)
        i <- i[is.finite(pos) & !is.na(q) & pos>=index$left[q] & pos<=index$right[q]]
      }
      hit <- if(length(i)) lines[i] else character()
    } else {
      snps <- sub('\t.*$','',lines)
      hit <- lines[snps %chin% sought]
    }
    if (length(hit)) {batches<-batches+1L; selected[[batches]]<-hit; nmatch<-nmatch+length(hit)}
    if (inspected-total>=5000000) {report(mode,'source data rows scanned',inspected,'window/SNP rows found',nmatch);total<-inspected}
  }
  report(mode, 'EOF REACHED. Source rows:',inspected,'extracted rows:',nmatch)
  if (inspected<1000000 || nmatch<5000) stop(mode,': suspiciously short/incomplete gzip or wrong index/window; do not use results')
  data <- fread(text=paste(c(hdr,unlist(selected,use.names=FALSE)),collapse='\n'),sep='\t',showProgress=FALSE)
  attr(data,'source_scanned')<-inspected
  data
}

htn_raw <- extract_stream(htn_path,'HTN')
htn <- htn_raw[,.(SNP=as.character(variant_id), chr_HTN=as.character(chromosome),
                  pos_HTN_GRCh37=as.numeric(base_pair_location), EA_HTN=toupper(as.character(effect_allele)),
                  OA_HTN=toupper(as.character(other_allele)), EAF_HTN=as.numeric(effect_allele_frequency),
                  beta_HTN=as.numeric(beta),se_HTN=as.numeric(standard_error),p_HTN=as.numeric(p_value))]
htn[,region:=index$region[match(chr_HTN,index$chr)]]
if (anyNA(htn$region)) stop('HTN region assignment missing')
save_tsv(htn,'02_HTN_raw_window_extracted.tsv')
report('HTN rows per audit region: ',paste(names(table(htn$region)),table(htn$region),collapse=' | '))
wanted <- unique(htn$SNP[grepl('^rs[0-9]+$',htn$SNP)])
wanted <- unique(c(wanted,as.character(index$index_snp),c('rs10509759','rs11191447','rs117159291','rs79780963')))
if(length(wanted)<9000) stop('Unexpectedly few HTN rsIDs in windows; inspect 02_HTN_raw_window_extracted.tsv')
mvp_raw <- extract_stream(mvp_path,'MVP',sought=wanted)
mvp <- mvp_raw[,.(SNP=as.character(SNP_ID),chr_MVP=as.character(chrom),pos_MVP=as.numeric(pos),
                  REF_MVP=toupper(as.character(ref)), ALT_MVP=toupper(as.character(alt)),
                  EA_MVP=toupper(as.character(ea)), EAF_MVP=as.numeric(af),
                  beta_MVP_raw=as.numeric(beta),se_MVP=as.numeric(sebeta),p_MVP=as.numeric(pval),
                  r2_MVP=as.numeric(r2),N_MVP=as.numeric(num_samples))]
save_tsv(mvp,'03_MVP_raw_SNP_extracted.tsv')
report('MVP rows per chromosome:',paste(names(table(mvp$chr_MVP)),table(mvp$chr_MVP),collapse=' | '))

# ID ambiguities: do not choose the row with smallest P; remove all copies of a duplicated rsID.
dup_htn <- unique(htn$SNP[duplicated(htn$SNP) | duplicated(htn$SNP,fromLast=TRUE)])
dup_mvp <- unique(mvp$SNP[duplicated(mvp$SNP) | duplicated(mvp$SNP,fromLast=TRUE)])
dup <- union(dup_htn,dup_mvp)
save_tsv(data.table(SNP=dup),'04_ambiguous_duplicated_rsIDs.tsv')
report('Duplicated identifiers: HTN',length(dup_htn),'MVP',length(dup_mvp),'union',length(dup))
a <- merge(htn[!SNP %chin% dup & grepl('^rs[0-9]+$',SNP)],
           mvp[!SNP %chin% dup & grepl('^rs[0-9]+$',SNP)],by='SNP')
a <- a[chr_HTN==chr_MVP]
if (nrow(a)<9000 || !setequal(unique(a$chr_HTN),c('10','18'))) stop('Unexpectedly few matched SNPs or missing region: check source files')
stopifnot(!anyDuplicated(a$SNP))
report('Matched, unique, same-chromosome SNPs:',nrow(a))

valid_base <- function(x) !is.na(x) & nchar(x)==1 & x %chin% c('A','C','G','T')
a[,snv_valid:=valid_base(EA_HTN)&valid_base(OA_HTN)&valid_base(EA_MVP)&
                 valid_base(REF_MVP)&valid_base(ALT_MVP)& EA_HTN!=OA_HTN & REF_MVP!=ALT_MVP &
                 (EA_MVP==REF_MVP | EA_MVP==ALT_MVP)]
a[,OA_MVP:=fifelse(EA_MVP==REF_MVP,ALT_MVP,
                         fifelse(EA_MVP==ALT_MVP,REF_MVP,NA_character_))]
a[,allele_relation:=fifelse(EA_HTN==EA_MVP & OA_HTN==OA_MVP,'same',
                            fifelse(EA_HTN==OA_MVP & OA_HTN==EA_MVP,'swapped','unresolved'))]
a[,palindromic:=paste0(EA_HTN,OA_HTN) %chin% c('AT','TA','CG','GC')]
a[,numeric_ok:=is.finite(beta_HTN)&is.finite(beta_MVP_raw)&is.finite(se_HTN)&is.finite(se_MVP)&
                se_HTN>0&se_MVP>0&is.finite(p_HTN)&is.finite(p_MVP)&
                p_HTN>=0&p_HTN<=1&p_MVP>=0&p_MVP<=1&is.finite(EAF_HTN)&
                is.finite(EAF_MVP)&EAF_HTN>=0&EAF_HTN<=1&EAF_MVP>=0&EAF_MVP<=1]
a[,beta_MVP_aligned:=fifelse(allele_relation=='same',beta_MVP_raw,
                            fifelse(allele_relation=='swapped',-beta_MVP_raw,NA_real_))]
a[,EAF_MVP_aligned:=fifelse(allele_relation=='same',EAF_MVP,
                            fifelse(allele_relation=='swapped',1-EAF_MVP,NA_real_))]
a[,freq_diff:=abs(EAF_HTN-EAF_MVP_aligned)]
a[,freq_ok:=is.finite(freq_diff) & freq_diff<=0.10]
a[,base_pass:=snv_valid & allele_relation %chin% c('same','swapped') & numeric_ok & freq_ok]
a[,regional_primary_pass:=base_pass & !palindromic]
a[,palindrome_sensitivity_pass:=base_pass & (!palindromic |
               (EAF_HTN<0.42 | EAF_HTN>0.58) & (EAF_MVP_aligned<0.42 | EAF_MVP_aligned>0.58))]
a[,Z_HTN:=beta_HTN/se_HTN]
a[,Z_MVP:=beta_MVP_aligned/se_MVP]
a[,direction_eligible:=p_HTN<1e-5 & beta_HTN!=0 & beta_MVP_aligned!=0]
a[,direction_concordant:=sign(beta_HTN)==sign(beta_MVP_aligned)]
save_tsv(a,'05_ALL_matched_SNP_QC_and_allele_orientation.tsv')

qc <- a[,.(n_matched=.N, n_biallelic_snv=sum(snv_valid),n_resolvable_alleles=sum(snv_valid&allele_relation!='unresolved'),
            n_numeric_valid=sum(numeric_ok),n_freq_compatible=sum(base_pass),
            n_palindromic_excluded=sum(base_pass&palindromic),
            n_primary_region=sum(regional_primary_pass),
            n_palindrome_sensitivity=sum(palindrome_sensitivity_pass)),by=region]
save_tsv(qc,'06_region_QC_counts.tsv')
report('REGIONAL QC:');report(paste(capture.output(print(qc)),collapse='\n'))

summarise <- function(subset,name) {
  x <- a[get(subset)==TRUE]
  if (!all(index$region %chin% x$region)) stop('No qualifying SNPs in one of the two regions: ',name)
  x[,.(analysis=name,n_regional=.N,pearson_Z=cor(Z_HTN,Z_MVP,method='pearson'),
       spearman_Z=cor(Z_HTN,Z_MVP,method='spearman'),
       n_direction_eligible=sum(direction_eligible),
       n_direction_concordant=sum(direction_eligible & direction_concordant),
       direction_fraction=if(sum(direction_eligible)) sum(direction_eligible & direction_concordant)/sum(direction_eligible) else NA_real_),by=region]
}
primary <- summarise('regional_primary_pass','raw_GWAS_nonpal_freqdiff_le_0.10')
sensitive <- summarise('palindrome_sensitivity_pass','raw_GWAS_freq_resolved_palindromes')
regional <- a[regional_primary_pass==TRUE]
save_tsv(regional,'07_region_QC_passed_primary.tsv')
save_tsv(rbindlist(list(primary,sensitive),use.names=TRUE), '08_main_and_palindrome_sensitivity_summary.tsv')
report('RAW GWAS PRIMARY REANALYSIS:');report(paste(capture.output(print(primary)),collapse='\n'))
report('PALINDROME SENSITIVITY:');report(paste(capture.output(print(sensitive)),collapse='\n'))

# Compare to prior archived *post hoc reanalysis*; discrepancies are REPORT ONLY, never used to modify QC.
expected <- data.table(region=c('chr10_region6','chr18_region3'),n_regional=c(5386L,5505L),
                       pearson_Z=c(0.641133589388029,0.2323640082918047),
                       spearman_Z=c(0.41436661665082336,0.20206429386660077),
                       n_direction_eligible=c(207L,70L),n_direction_concordant=c(207L,70L))
cmp <- merge(expected,primary,by='region',suffixes=c('_previous','_raw_R'),all=TRUE)
cmp[,counts_match:=n_regional_previous==n_regional_raw_R &
                 n_direction_eligible_previous==n_direction_eligible_raw_R &
                 n_direction_concordant_previous==n_direction_concordant_raw_R]
cmp[,corr_match:=abs(pearson_Z_previous-pearson_Z_raw_R)<1e-6 &
               abs(spearman_Z_previous-spearman_Z_raw_R)<1e-6]
save_tsv(cmp,'09_compare_previous_reanalysis_VS_raw_GWAS_R.tsv')
report('COMPARISON WITH PREVIOUS REANALYSIS:');report(paste(capture.output(print(cmp)),collapse='\n'))
if (!all(cmp$counts_match %in% TRUE) || !all(cmp$corr_match %in% TRUE)) {
  report('REVIEW REQUIRED: raw R does not match previous reanalysis; DO NOT revise Figure 3/manuscript yet.')
} else report('PASS: raw GWAS R reproduces the previously archived post hoc regional summary to defined tolerance.')

# Representative variants (includes palindromic index SNPs if allele and frequency checks permit).
rep_snps <- c('rs1060240','rs10509759','rs11191447','rs1788799')
rep <- a[SNP %chin% rep_snps]
save_tsv(rep,'10_four_representative_variants_raw_GWAS.tsv')
if (!setequal(rep$SNP,rep_snps)) stop('Representative SNP missing after source merge: inspect duplicates/coordinates')
if (!all(rep$base_pass)) report('REVIEW REQUIRED: one or more representative variants failed allele/frequency QC.')
lead <- fread(lead_path)
lead_check <- merge(rep[,.(SNP,region,EA_HTN,OA_HTN,EA_MVP,palindromic,beta_HTN,p_HTN,
                           beta_MVP_aligned,se_MVP,p_MVP,EAF_HTN,EAF_MVP_aligned,freq_diff,base_pass)],
                    lead[,.(SNP=topSNP,beta_HTN_previous=HTN_beta,p_HTN_previous=HTN_p,
                            beta_MVP_previous=MVP_beta_harmonized,p_MVP_previous=MVP_p)],by='SNP',all=TRUE)
lead_check[,p_HTN_log10_diff:=abs(log10(pmax(p_HTN, .Machine$double.xmin))-
                                   log10(pmax(p_HTN_previous,.Machine$double.xmin)))]
lead_check[,p_MVP_log10_diff:=abs(log10(pmax(p_MVP,.Machine$double.xmin))-
                                   log10(pmax(p_MVP_previous,.Machine$double.xmin)))]
lead_check[,lead_effect_match:=abs(beta_HTN-beta_HTN_previous)<1e-6 &
              abs(beta_MVP_aligned-beta_MVP_previous)<1e-6 &
              p_HTN_log10_diff<1e-3 & p_MVP_log10_diff<1e-3]
save_tsv(lead_check,'11_representative_variants_comparison.tsv')
report('REPRESENTATIVE SNP COMPARISON:');report(paste(capture.output(print(lead_check)),collapse='\n'))
if (!all(lead_check$lead_effect_match %in% TRUE)) report('REVIEW REQUIRED: representative SNP disagreement with archived table.')

# Exploration only: strongest MVP association among SAME eligible ±1Mb primary regional SNPs.
minima <- regional[order(p_MVP,SNP),.SD[1],by=region][,.(region,SNP,pos_HTN_GRCh37,p_MVP,p_HTN,
                       beta_MVP_aligned,beta_HTN,EA_HTN,EA_MVP)]
minima[,distance_from_index_bp:=abs(pos_HTN_GRCh37-index$center[match(region,index$region)])]
save_tsv(minima,'12_exploratory_regional_minima_NOT_primary_tests.tsv')
report('EXPLORATORY MINIMA:');report(paste(capture.output(print(minima)),collapse='\n'))
chr18_detail <- a[SNP %chin% c('rs1788799','rs117159291')]
save_tsv(chr18_detail,'13_chr18_index_vs_remote_peak_no_LD_inference.tsv')

hist <- fread(hist_path)
save_tsv(hist,'14_historical_summary_READ_ONLY.tsv')
manifest <- data.table(input=names(inputs),path=unname(inputs),size_bytes=as.numeric(file.info(inputs)$size),
                       md5=unname(tools::md5sum(inputs)))
save_tsv(manifest,'15_input_provenance_size_MD5.tsv')
writeLines(capture.output(sessionInfo()),file.path(out,'16_R_sessionInfo.txt'))
writeLines(c('Figure 3 raw-GWAS final audit, POST HOC; no GWAS rerun.',
             'Region: GRCh37 index SNP +/-1Mb; match by rsID+same chromosome+allele set; exclude duplicate IDs,',
             'exclude palindrome in primary regional set, <=0.10 matched effect-allele-frequency difference.',
             'Direction subset: HTN p<1e-5, both effects !=0, denominator recorded separately from total regional SNPs.',
             'Palindrome sensitivity requires both aligned allele frequencies outside [0.42,0.58].',
             'Historical 2026-09-13 regional protocol is NOT thereby recovered.',
             'Do not use regional-minimum P as a 2-primary-variant Bonferroni test.',
             'No rs1788799-rs117159291 LD calculation was performed in this script.',
             'Review 09 and 11 comparison files before updating figure/manuscript.'),
           file.path(out,'README_methods_and_limits.txt'))
report('FINISHED. Review 08, 09, 11, 12, 13 before any manuscript edits. Path:',out)
