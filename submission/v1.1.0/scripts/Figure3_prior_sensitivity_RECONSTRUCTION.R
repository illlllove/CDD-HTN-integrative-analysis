# ============================================================
# Figure 3
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

data_dir <- file.path(root, "frozen_data", "Figure3")
out_dir  <- file.path(root, "figures", "reconstructed")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

count <- fread(file.path(
  data_dir,
  "Figure3A_candidate_count_by_prior.tsv"
))

post <- fread(file.path(
  data_dir,
  "Figure3BC_HTN_CDD_H3_H4_by_prior.tsv"
))

minh4 <- fread(file.path(
  data_dir,
  "Figure3D_four_candidates_minimum_H4_by_prior.tsv"
))

prior_levels <- c("1e-06", "1e-05", "5e-05")
prior_labels <- c(
  "1 x 10^-6",
  "1 x 10^-5",
  "5 x 10^-5"
)

count[, p12_f := factor(
  as.character(p12),
  levels = prior_levels,
  labels = prior_labels
)]

post[, p12_f := factor(
  as.character(p12),
  levels = prior_levels,
  labels = prior_labels
)]

minh4[, p12_f := factor(
  as.character(p12),
  levels = prior_levels,
  labels = prior_labels
)]

# ------------------------------------------------------------
# Panel A
# ------------------------------------------------------------

pA <- ggplot(count, aes(p12_f, n_candidates)) +
  geom_col(width = 0.65) +
  geom_text(
    aes(label = n_candidates),
    vjust = -0.4,
    size = 3.5
  ) +
  labs(
    title = "A  Candidate retention by shared-variant prior",
    subtitle = "412 exploratory probes; all three pairwise H4 >= 0.50",
    x = "p12",
    y = "Retained candidates"
  ) +
  theme_bw(base_size = 9)

# ------------------------------------------------------------
# Panels B-C
# ------------------------------------------------------------

make_posterior_panel <- function(gene, panel_letter) {

  d <- post[Gene == gene]

  dl <- rbind(
    d[, .(p12_f, Hypothesis = "H3", Posterior = H3)],
    d[, .(p12_f, Hypothesis = "H4", Posterior = H4)]
  )

  ggplot(
    dl,
    aes(
      p12_f,
      Posterior,
      group = Hypothesis,
      linetype = Hypothesis,
      shape = Hypothesis
    )
  ) +
    annotate(
      "rect",
      xmin = 1.75,
      xmax = 2.25,
      ymin = -Inf,
      ymax = Inf,
      alpha = 0.08
    ) +
    geom_line() +
    geom_point(size = 2.4) +
    geom_text(
      aes(label = sprintf("%.3f", Posterior)),
      vjust = -0.65,
      size = 3
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      title = paste0(
        panel_letter,
        "  ",
        gene,
        ": HTN-CDD posterior sensitivity"
      ),
      subtitle = "Grey band indicates the primary prior",
      x = "p12",
      y = "Posterior probability"
    ) +
    theme_bw(base_size = 9)
}

pB <- make_posterior_panel("NT5C2", "B")
pC <- make_posterior_panel("RMC1", "C")

# ------------------------------------------------------------
# Panel D
# ------------------------------------------------------------

gene_order <- c(
  "NT5C2",
  "RMC1",
  "CNNM2",
  "MARCKSL1P1"
)

minh4[, Gene := factor(Gene, levels = rev(gene_order))]

pD <- ggplot(
  minh4,
  aes(
    p12_f,
    Gene,
    fill = min_H4
  )
) +
  geom_tile() +
  geom_text(
    aes(label = sprintf("%.3f", min_H4)),
    size = 3
  ) +
  scale_fill_gradient(
    low = "white",
    high = "black",
    limits = c(0, 1)
  ) +
  labs(
    title = "D  Minimum pairwise H4 by prior",
    subtitle = "Weakest H4 across eQTL-HTN, eQTL-CDD and HTN-CDD",
    x = "p12",
    y = NULL,
    fill = "Minimum H4"
  ) +
  theme_bw(base_size = 9)

final_plot <- (pA | pD) / (pB | pC)

ggsave(
  file.path(
    out_dir,
    "Figure3_prior_sensitivity_RECONSTRUCTED.pdf"
  ),
  final_plot,
  width = 11,
  height = 8
)

message("Figure 3 reconstruction complete.")
