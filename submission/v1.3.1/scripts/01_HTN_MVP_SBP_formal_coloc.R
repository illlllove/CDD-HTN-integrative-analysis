suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

# Usage:
# Rscript 01_HTN_MVP_SBP_formal_coloc.R <07_region_QC_passed_primary.tsv> <output_dir>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop('Usage: Rscript 01_HTN_MVP_SBP_formal_coloc.R <input.tsv> <output_dir>')
input_file <- args[1]
out_dir <- args[2]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dat <- fread(input_file)
req <- c('SNP','region','pos_HTN_GRCh37','beta_HTN','se_HTN','beta_MVP_aligned','se_MVP')
miss <- setdiff(req, names(dat))
if (length(miss)) stop('Missing columns: ', paste(miss, collapse=', '))

dat <- dat[!is.na(SNP) & !is.na(region) & is.finite(beta_HTN) & is.finite(se_HTN) & se_HTN>0 &
             is.finite(beta_MVP_aligned) & is.finite(se_MVP) & se_MVP>0]

run_one <- function(x, reg, p12) {
  d1 <- list(beta=as.numeric(x$beta_HTN), varbeta=as.numeric(x$se_HTN)^2,
             snp=as.character(x$SNP), position=as.numeric(x$pos_HTN_GRCh37), type='cc')
  d2 <- list(beta=as.numeric(x$beta_MVP_aligned), varbeta=as.numeric(x$se_MVP)^2,
             snp=as.character(x$SNP), position=as.numeric(x$pos_HTN_GRCh37), type='quant', sdY=1)
  fit <- coloc.abf(d1, d2, p1=1e-4, p2=1e-4, p12=p12)
  s <- fit$summary
  data.table(region=reg,p12=p12,nsnps=as.integer(s['nsnps']),
             PP.H0=as.numeric(s['PP.H0.abf']),PP.H1=as.numeric(s['PP.H1.abf']),
             PP.H2=as.numeric(s['PP.H2.abf']),PP.H3=as.numeric(s['PP.H3.abf']),
             PP.H4=as.numeric(s['PP.H4.abf']))
}

regions <- c('chr10_region6','chr18_region3')
priors <- c(1e-6,1e-5,5e-5)
res <- rbindlist(lapply(regions, function(reg) {
  x <- dat[region==reg]
  rbindlist(lapply(priors, function(p12) run_one(x, reg, p12)))
}))
res[, `:=`(H4_minus_H3=PP.H4-PP.H3,
           H4_over_H3=fifelse(PP.H3>0,PP.H4/PP.H3,Inf),
           H4_ge_0.50=PP.H4>=0.50)]
fwrite(res, file.path(out_dir,'HTN_MVP_SBP_formal_coloc_ALL_PRIORS.tsv'), sep='\t')
capture.output(sessionInfo(), file=file.path(out_dir,'sessionInfo.txt'))
print(res)
