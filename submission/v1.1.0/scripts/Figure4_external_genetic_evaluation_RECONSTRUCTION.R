# ============================================================
# Figure 4
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

data_dir <- file.path(root, "frozen_data", "Figure4")
out_dir  <- file.path(root, "figures", "reconstructed")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

corr <- fread(file.path(
  data_dir,
  "Figure4A_MVP_regional_correlation.tsv"
))

dirc <- fread(file.path(
  data_dir,
  "Figure4B_MVP_direction_concordance.tsv"
))

h4 <- fread(file.path(
  data_dir,
  "Figure4C_discovery_vs_UKB_M50_H4.tsv"
))

# ------------------------------------------------------------
# Panel A
# ------------------------------------------------------------

pA <- ggplot(
  corr,
  aes(
    Metric,
    Correlation,
    fill = Region
  )
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65
  ) +
  geom_text(
    aes(label = sprintf("%.3f", Correlation)),
    position = position_dodge(width = 0.75),
    vjust = -0.4,
    size = 3
  ) +
  coord_cartesian(ylim = c(0, 0.72)) +
  labs(
    title = "A  Regional concordance with MVP SBP",
    subtitle = "Aligned Z-score correlations in conservative reconstructed SNP sets",
    x = NULL,
    y = "Correlation"
  ) +
  theme_bw(base_size = 9)

# ------------------------------------------------------------
# Panel B
# ------------------------------------------------------------

dirc[, percent := proportion * 100]

pB <- ggplot(
  dirc,
  aes(
    Region,
    percent
  )
) +
  geom_col(width = 0.6) +
  geom_text(
    aes(label = label),
    vjust = -0.5,
    size = 3.5
  ) +
  coord_cartesian(ylim = c(0, 107)) +
  labs(
    title = "B  Direction concordance with MVP SBP",
    subtitle = "Discovery HTN P < 1e-5; SNPs are LD-correlated",
    x = NULL,
    y = "Concordant direction (%)"
  ) +
  theme_bw(base_size = 9)

# ------------------------------------------------------------
# Panel C
# ------------------------------------------------------------

gene_order <- c(
  "NT5C2",
  "RMC1",
  "CNNM2",
  "MARCKSL1P1"
)

h4[, Gene := factor(Gene, levels = rev(gene_order))]

h4_long <- rbind(
  h4[, .(
    Gene,
    Dataset = "Discovery CDD",
    H4 = H4_discovery
  )],
  h4[, .(
    Gene,
    Dataset = "UKB M50",
    H4 = H4_UKB
  )]
)

pC <- ggplot(
  h4_long,
  aes(
    H4,
    Gene,
    group = Gene
  )
) +
  geom_vline(
    xintercept = 0.50,
    linetype = 2
  ) +
  geom_line(
    aes(group = Gene)
  ) +
  geom_point(
    aes(shape = Dataset),
    size = 2.7
  ) +
  coord_cartesian(xlim = c(0, 0.82)) +
  labs(
    title = "C  CDD-side colocalization in discovery and UKB M50",
    subtitle = "Primary prior p12 = 1e-5; dashed line marks H4 = 0.50",
    x = "PP.H4",
    y = NULL
  ) +
  theme_bw(base_size = 9)

final_plot <- (pA | pB) / pC

ggsave(
  file.path(
    out_dir,
    "Figure4_external_genetic_evaluation_RECONSTRUCTED.pdf"
  ),
  final_plot,
  width = 10,
  height = 7
)

message("Figure 4 reconstruction complete.")
