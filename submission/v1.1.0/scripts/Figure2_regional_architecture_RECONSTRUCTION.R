# ============================================================
# Figure 2
# Submission-specific reconstruction from frozen derived tables
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)

if (length(script_arg) == 1) {
  script_path <- normalizePath(sub("^--file=", "", script_arg),
                               winslash = "/", mustWork = FALSE)
  root <- dirname(dirname(script_path))
} else {
  root <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

data_dir <- file.path(root, "frozen_data", "Figure2")
out_dir  <- file.path(root, "figures", "reconstructed")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

chr10 <- fread(file.path(
  data_dir,
  "Figure2_chr10_NT5C2_regional_master.tsv"
))

chr18 <- fread(file.path(
  data_dir,
  "Figure2_chr18_NPC1_RMC1_regional_master.tsv"
))

genes <- fread(file.path(
  data_dir,
  "Figure2_GRCh37_gene_ranges_FINAL.tsv"
))

# ------------------------------------------------------------
# helper
# ------------------------------------------------------------

regional_panel <- function(dat, chr, title,
                           global_h4,
                           cs_n,
                           gene_dat) {

  dat[, pos_mb := position / 1e6]

  p1 <- ggplot(dat, aes(pos_mb, minuslog10P_HTN)) +
    geom_point(size = 0.8, alpha = 0.75) +
    labs(
      title = title,
      y = "-log10(P) HTN",
      x = NULL
    ) +
    theme_bw(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )

  p2 <- ggplot(dat, aes(pos_mb, minuslog10P_CDD)) +
    geom_point(size = 0.8, alpha = 0.75) +
    labs(
      y = "-log10(P) CDD",
      x = NULL
    ) +
    theme_bw(base_size = 9) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )

  p3 <- ggplot(dat, aes(pos_mb, SNP.PP.H4)) +
    geom_point(size = 0.9, alpha = 0.8) +
    labs(
      subtitle = paste0(
        "Global PP.H4 = ", global_h4,
        "; 95% conditional H4 CS = ", cs_n, " SNPs"
      ),
      y = "Conditional SNP.PP.H4",
      x = NULL
    ) +
    theme_bw(base_size = 9) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )

  gd <- gene_dat[chromosome == chr]
  gd[, start_mb := start / 1e6]
  gd[, end_mb   := end / 1e6]
  gd[, y := seq_len(.N)]

  p4 <- ggplot(gd) +
    geom_segment(
      aes(
        x = start_mb,
        xend = end_mb,
        y = y,
        yend = y
      ),
      linewidth = 1
    ) +
    geom_text(
      aes(
        x = (start_mb + end_mb) / 2,
        y = y + 0.2,
        label = gene
      ),
      size = 3
    ) +
    scale_y_continuous(NULL, breaks = NULL) +
    labs(
      x = "Genomic position (Mb, GRCh37)"
    ) +
    theme_bw(base_size = 9)

  p1 / p2 / p3 / p4 +
    plot_layout(heights = c(1, 1, 1, 0.7))
}


p10 <- regional_panel(
  chr10,
  chr = 10,
  title = "A  Chromosome 10 candidate region",
  global_h4 = "0.756",
  cs_n = 39,
  gene_dat = genes
)

p18 <- regional_panel(
  chr18,
  chr = 18,
  title = "B  Chromosome 18 NPC1-RMC1 region",
  global_h4 = "0.647",
  cs_n = 42,
  gene_dat = genes
)

final_plot <- p10 | p18

ggsave(
  file.path(
    out_dir,
    "Figure2_regional_architecture_RECONSTRUCTED.pdf"
  ),
  final_plot,
  width = 12,
  height = 8.5
)

message("Figure 2 reconstruction complete.")
