suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

# This script expects the eight frozen formal candidate input tables used for the
# eQTL-HTN and eQTL-CDD comparisons. Set ROOT to the local project path.
ROOT <- Sys.getenv('CDD_HTN_ROOT', unset='D:/Bioinformatics/CDD_HTN')
OUT <- file.path(ROOT,'10_sdY_sensitivity')
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

inputs <- data.table(
  Gene=c('NT5C2','NT5C2','CNNM2','CNNM2','MARCKSL1P1','MARCKSL1P1','RMC1','RMC1'),
  Pair=c('eQTL-HTN','eQTL-CDD','eQTL-HTN','eQTL-CDD','eQTL-HTN','eQTL-CDD','eQTL-HTN','eQTL-CDD'),
  File=c(
    file.path(ROOT,'09_manuscript_audit/NT5C2_FORMAL_RECONSTRUCTION/02_NT5C2_EQTL_HTN_SNP_input.tsv'),
    file.path(ROOT,'09_manuscript_audit/NT5C2_FORMAL_RECONSTRUCTION/03_NT5C2_EQTL_CDD_SNP_input.tsv'),
    file.path(ROOT,'09_manuscript_audit/FORMAL_COLOC_SNP_EXPORT/CNNM2_EQTL_HTN_SNP_input.tsv'),
    file.path(ROOT,'09_manuscript_audit/FORMAL_COLOC_SNP_EXPORT/CNNM2_EQTL_CDD_SNP_input.tsv'),
    file.path(ROOT,'09_manuscript_audit/FORMAL_COLOC_SNP_EXPORT/MARCKSL1P1_EQTL_HTN_SNP_input.tsv'),
    file.path(ROOT,'09_manuscript_audit/FORMAL_COLOC_SNP_EXPORT/MARCKSL1P1_EQTL_CDD_SNP_input.tsv'),
    file.path(ROOT,'09_manuscript_audit/FORMAL_COLOC_SNP_EXPORT/RMC1_EQTL_HTN_SNP_input.tsv'),
    file.path(ROOT,'09_manuscript_audit/FORMAL_COLOC_SNP_EXPORT/RMC1_EQTL_CDD_SNP_input.tsv')
  )
)

# Column mappings differ for NT5C2 vs the other archived exports.
read_input <- function(gene, path) {
  d <- fread(path)
  if (gene=='NT5C2') {
    list(snp=d$SNP, eb=d$eqtl_b, ese=d$eqtl_se, gb=d$GWAS_b_aligned, gse=d$GWAS_se,
         keep=!d$palindromic & !d$mismatch & !d$freq_fail)
  } else {
    list(snp=d$SNP, eb=d$b, ese=d$SE, gb=d$GW_b_aligned, gse=d$GW_se,
         keep=!d$palindromic & !d$mismatch & !d$freq_fail)
  }
}

run_pair <- function(gene,pair,path,sdY,p12=1e-5) {
  x <- read_input(gene,path); k <- x$keep
  d1 <- list(beta=x$eb[k], varbeta=x$ese[k]^2, snp=as.character(x$snp[k]), type='quant', sdY=sdY)
  d2 <- list(beta=x$gb[k], varbeta=x$gse[k]^2, snp=as.character(x$snp[k]), type='cc')
  fit <- coloc.abf(d1,d2,p1=1e-4,p2=1e-4,p12=p12)
  data.table(Gene=gene,Pair=pair,sdY=sdY,n=as.integer(fit$summary['nsnps']),H4=as.numeric(fit$summary['PP.H4.abf']))
}

res <- rbindlist(lapply(seq_len(nrow(inputs)), function(i) {
  rbindlist(lapply(c(0.8,1.0,1.2), function(s) run_pair(inputs$Gene[i],inputs$Pair[i],inputs$File[i],s)))
}))
fwrite(res,file.path(OUT,'sdY_sensitivity_REFERENCE_p12_1e-5_long.tsv'),sep='\t')
print(res)
