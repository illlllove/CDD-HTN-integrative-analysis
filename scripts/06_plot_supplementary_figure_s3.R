#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(scales)
  library(ragg)
})

root <- Sys.getenv("CDD_HTN_ROOT", unset = ".")
source_dir <- file.path(
  root,
  "07_integration",
  "S3_bulk_plotData_extraction_FINAL_2026-09-19_v3"
)
validation_file <- file.path(
  root,
  "05_expression",
  "results",
  "GSE34095_frozen_candidate_secondary_validation_2026-09-13.tsv"
)
out_dir <- Sys.getenv(
  "CDD_HTN_OUTPUT_DIR",
  unset = file.path(root, "outputs", "figures")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

candidate_levels <- c("NT5C2", "CNNM2", "MARCKSL1P1", "RMC1")
dataset_levels <- c("GSE153761", "GSE28360", "GSE34095")
dataset_labels <- c(
  GSE153761 = "Primary CDD\nGSE153761",
  GSE28360 = "Hypertension\nGSE28360",
  GSE34095 = "Disc-related context\nGSE34095"
)

coverage <- read_tsv(
  file.path(source_dir, "S3_candidate_dataset_coverage_matrix.tsv"),
  show_col_types = FALSE
) %>%
  mutate(
    gene = factor(gene, levels = rev(candidate_levels)),
    dataset = factor(dataset, levels = dataset_levels),
    tile_text = if_else(status == "Evaluable", "Eval.", "NC"),
    tile_class = if_else(status == "Evaluable", "Evaluable", "Not covered")
  )

stats <- read_tsv(
  file.path(source_dir, "S3_mapping_driven_candidate_statistics_ALL.tsv"),
  show_col_types = FALSE,
  na = c("", "NA", '""')
) %>%
  mutate(
    gene = factor(gene, levels = rev(candidate_levels)),
    dataset = factor(dataset, levels = dataset_levels),
    dataset_display = factor(dataset_labels[as.character(dataset)],
                             levels = unname(dataset_labels)),
    nominal_class = if_else(stored_P < 0.05, "Nominal P < 0.05", "Nominal P >= 0.05"),
    p_label = paste0("P=", format.pval(stored_P, digits = 2, eps = 0.001)),
    label_x = if_else(stored_logFC >= 0, stored_logFC - 0.10, stored_logFC + 0.10),
    label_vjust = -0.95,
    label_x = case_when(
      dataset == "GSE34095" & feature_id == "209155_s_at" ~ 0.13,
      dataset == "GSE34095" & feature_id == "206818_s_at" ~ 0.08,
      dataset == "GSE34095" & feature_id == "209874_x_at" ~ -0.24,
      dataset == "GSE34095" & feature_id == "221190_s_at" ~ 0.07,
      TRUE ~ label_x
    ),
    label_vjust = case_when(
      dataset == "GSE34095" & feature_id == "209874_x_at" ~ 1.35,
      TRUE ~ label_vjust
    )
  )

validation <- read_tsv(validation_file, show_col_types = FALSE)
nt5c2 <- validation %>% filter(candidate_gene == "NT5C2", probe_id == "209155_s_at")
stopifnot(nrow(nt5c2) == 1)
stopifnot(abs(nt5c2$P.Value - 0.025363300206407) < 1e-12)
stopifnot(abs(nt5c2$P_BH_candidate_probe - 0.101453200825628) < 1e-12)
stopifnot(abs(nt5c2$BH_genomewide - 0.917002857797789) < 1e-12)

context <- crossing(
  gene = factor(candidate_levels, levels = rev(candidate_levels)),
  context = factor(c("Primary CDD\nGSE153761", "Hypertension\nGSE28360",
                     "Disc-related context\nGSE34095", "Overall\nbulk context"),
                   levels = c("Primary CDD\nGSE153761", "Hypertension\nGSE28360",
                              "Disc-related context\nGSE34095", "Overall\nbulk context"))
) %>%
  mutate(label = "", class = "No support")

set_context <- function(gene_name, context_name, label_value, class_value) {
  idx <- context$gene == gene_name & context$context == context_name
  context$label[idx] <<- label_value
  context$class[idx] <<- class_value
}

set_context("NT5C2", "Primary CDD\nGSE153761", "No DE\nsupport", "No support")
set_context("NT5C2", "Hypertension\nGSE28360", "No DE\nsupport", "No support")
set_context("NT5C2", "Disc-related context\nGSE34095",
            "Nominal only\nnot BH-significant\nage-confounded", "Nominal only")
set_context("NT5C2", "Overall\nbulk context", "Contextual only\nnot confirmatory", "Contextual")

set_context("CNNM2", "Primary CDD\nGSE153761", "No DE\nsupport", "No support")
set_context("CNNM2", "Hypertension\nGSE28360", "No DE\nsupport", "No support")
set_context("CNNM2", "Disc-related context\nGSE34095", "No DE\nsupport", "No support")
set_context("CNNM2", "Overall\nbulk context", "No bulk\nsupport", "Overall no support")

set_context("MARCKSL1P1", "Primary CDD\nGSE153761", "Not\ncovered", "Not covered")
set_context("MARCKSL1P1", "Hypertension\nGSE28360", "Not\ncovered", "Not covered")
set_context("MARCKSL1P1", "Disc-related context\nGSE34095", "Not\ncovered", "Not covered")
set_context("MARCKSL1P1", "Overall\nbulk context", "Not\nevaluable", "Not evaluable")

set_context("RMC1", "Primary CDD\nGSE153761", "Not\ncovered", "Not covered")
set_context("RMC1", "Hypertension\nGSE28360", "No DE\nsupport", "No support")
set_context("RMC1", "Disc-related context\nGSE34095", "No DE\nsupport", "No support")
set_context("RMC1", "Overall\nbulk context", "No bulk\nsupport", "Overall no support")

base_theme <- theme_minimal(base_size = 8.5) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(colour = "#222222"),
    axis.text = element_text(colour = "#222222"),
    strip.text = element_text(face = "bold", size = 8.2),
    plot.tag = element_text(face = "bold", size = 15),
    plot.tag.position = c(0, 1),
    plot.margin = margin(7, 7, 5, 7)
  )

p_a <- ggplot(coverage, aes(dataset, gene, fill = tile_class)) +
  geom_tile(colour = "white", linewidth = 1.1) +
  geom_text(aes(label = tile_text, colour = tile_class), fontface = "bold", size = 3.0) +
  scale_fill_manual(values = c("Evaluable" = "#087FB9", "Not covered" = "#DEDEDE"), guide = "none") +
  scale_colour_manual(values = c("Evaluable" = "white", "Not covered" = "#222222"), guide = "none") +
  scale_x_discrete(labels = c("153761", "28360", "34095")) +
  labs(x = "GSE accession", y = NULL) +
  coord_equal() +
  base_theme +
  theme(axis.text.x = element_text(size = 7.8), axis.text.y = element_text(face = "italic"))

p_b <- ggplot(stats, aes(stored_logFC, gene)) +
  geom_vline(xintercept = 0, linewidth = 0.45, colour = "#777777") +
  geom_point(aes(fill = nominal_class), shape = 21, size = 2.8, stroke = 0.7, colour = "black") +
  geom_text(aes(x = label_x, label = p_label, vjust = label_vjust),
            size = 2.35, hjust = 0.5) +
  facet_wrap(~dataset_display, nrow = 1) +
  scale_fill_manual(values = c("Nominal P < 0.05" = "#E97619", "Nominal P >= 0.05" = "white"),
                    name = NULL) +
  scale_x_continuous(limits = c(-0.34, 0.38), breaks = c(-0.2, 0, 0.2)) +
  labs(
    x = "Stored logFC",
    y = NULL,
    caption = paste0(
      "GSE34095 NT5C2: nominal P=0.025; BH candidate-probe P=0.101;\n",
      "BH genome-wide P=0.917 (contextual only; not confirmatory)"
    )
  ) +
  base_theme +
  theme(
    panel.grid.major.x = element_line(colour = "#E1E1E1", linewidth = 0.35),
    axis.text.y = element_text(face = "italic"),
    legend.position = "bottom",
    legend.key.height = unit(3.2, "mm"),
    plot.caption = element_text(size = 6.9, hjust = 0.5, colour = "#333333",
                                lineheight = 1.0, margin = margin(t = 3))
  )

p_c <- ggplot(context, aes(context, gene, fill = class)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = label), size = 2.55, lineheight = 0.92) +
  scale_fill_manual(
    values = c(
      "No support" = "#D8D8D8",
      "Nominal only" = "#F0A202",
      "Contextual" = "#F6E83B",
      "Not covered" = "#EEEEEE",
      "Not evaluable" = "#FFFFFF",
      "Overall no support" = "#BDBDBD"
    ),
    guide = "none"
  ) +
  labs(x = NULL, y = NULL) +
  base_theme +
  theme(
    axis.text.x = element_text(size = 7.7),
    axis.text.y = element_text(face = "italic"),
    plot.margin = margin(4, 7, 7, 7)
  )

final_plot <- ((p_a + p_b) + plot_layout(widths = c(0.75, 1.75))) /
  p_c +
  plot_layout(heights = c(1.05, 1.0)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 15))

pdf_file <- file.path(out_dir, "Supplementary_Figure_S3_FINAL.pdf")
tiff_file <- file.path(out_dir, "Supplementary_Figure_S3_FINAL_600dpi.tiff")
png_file <- file.path(out_dir, "Supplementary_Figure_S3_FINAL_preview.png")

ggsave(pdf_file, final_plot, width = 180, height = 140, units = "mm",
       device = cairo_pdf, bg = "white")
ragg::agg_tiff(tiff_file, width = 180, height = 140, units = "mm",
               res = 600, compression = "lzw", background = "white")
print(final_plot)
dev.off()
ragg::agg_png(png_file, width = 180, height = 140, units = "mm",
              res = 180, background = "white")
print(final_plot)
dev.off()

message("Created:")
message(pdf_file)
message(tiff_file)
message(png_file)
