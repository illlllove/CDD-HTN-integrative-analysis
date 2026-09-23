# CDD–HTN exploratory shared-SMR screen and saved-coloc selection reconstruction
# Base R only; NO GWAS, SMR or colocalization computation is performed.
# From the extracted package root:
# Rscript 06_code/reproduce_SMR_and_coloc.R
# Optionally pass package root as the sole argument.
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1]] else "."
root <- normalizePath(root, mustWork = TRUE)
out <- file.path(root, "05_results_recalculated", "R_reconstruction")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cdd_path <- file.path(root,"01_primary_inputs","CDD_eQTLGen_all.smr")
htn_path <- file.path(root,"01_primary_inputs","HTN_eQTLGen_all.smr")
batch_path <- file.path(root,"03_historical_selection","batch_coloc_412_all_priors.tsv")
freeze_path <- file.path(root,"03_historical_selection","candidate_freeze_2026-09-13.tsv")
read_smr <- function(path) {
  x <- read.table(path, header=TRUE, sep="", quote="", comment.char="",
                  check.names=FALSE, stringsAsFactors=FALSE, na.strings=c("NA",""))
  stopifnot(!anyDuplicated(x$probeID), all(!is.na(x$p_SMR)),
            all(is.finite(x$p_SMR)), all(x$p_SMR >= 0 & x$p_SMR <= 1))
  x
}
cdd <- read_smr(cdd_path); htn <- read_smr(htn_path)
common <- merge(cdd, htn, by="probeID", suffixes=c("_CDD","_HTN"), sort=TRUE)
stopifnot(nrow(common)>0)
common$Pshared <- pmax(common$p_SMR_CDD, common$p_SMR_HTN)
common$q_BH_shared <- p.adjust(common$Pshared, method="BH")
common$nominal_both <- common$Pshared < 0.05
common$HEIDI_both_gt_001 <- !is.na(common$p_HEIDI_CDD) &
   !is.na(common$p_HEIDI_HTN) &
   common$p_HEIDI_CDD > 0.01 & common$p_HEIDI_HTN > 0.01
common <- common[order(common$Pshared,common$probeID),]
common$p_shared_order <- seq_len(nrow(common))
write.table(common, file.path(out,"shared_SMR_all_R.tsv"), sep="\t", quote=FALSE, row.names=FALSE, na="NA")
nominal <- common[common$nominal_both,]
write.table(nominal,file.path(out,"shared_SMR_nominal_R.tsv"),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
batch <- read.delim(batch_path, check.names=FALSE, stringsAsFactors=FALSE)
stopifnot(!anyDuplicated(paste(batch$probeID, batch$pair, batch$p12)))
expected_pairs <- c("EQTL_HTN","EQTL_CDD","HTN_CDD")
select_at_prior <- function(prior) {
  sub <- batch[abs(batch$p12-prior)<1e-12,]
  good <- sub[!is.na(sub$H4) & sub$H4 >= 0.50 & sub$pair %in% expected_pairs,]
  tab <- table(good$probeID)
  sort(names(tab)[tab == 3L])
}
prior <- c(1e-6,1e-5,5e-5)
selected <- lapply(prior, select_at_prior)
freeze <- read.delim(freeze_path,check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(identical(selected[[2]],sort(as.character(freeze$probeID))))
counts <- data.frame(p12=format(prior,scientific=TRUE),
                     number_meeting_all_three_H4_ge_0.50=lengths(selected))
write.table(counts,file.path(out,"coloc_prior_sensitivity_counts_R.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
cat("PASS: CDD",nrow(cdd),"HTN",nrow(htn),"joint",nrow(common),
    "nominal",nrow(nominal),"BH_significant",sum(common$q_BH_shared<0.05),"\n")
cat("PASS: prior counts",paste(lengths(selected),collapse=","),
    "freeze-primary-match",identical(selected[[2]],sort(as.character(freeze$probeID))),"\n")
cat("IMPORTANT: Cannot prove original threshold prespecification or historical analysis chronology.\n")
