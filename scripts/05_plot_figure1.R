#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(grid)
  library(ragg)
})

root <- Sys.getenv("CDD_HTN_ROOT", unset = ".")
out_dir <- Sys.getenv(
  "CDD_HTN_OUTPUT_DIR",
  unset = file.path(root, "outputs", "figures")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

blue_fill <- "#DCEAF5"
blue_border <- "#6E8798"
brown_fill <- "#F4EEE9"
brown_border <- "#9A806B"
grey_fill <- "#F7F7F7"
dark <- "#222222"
mid <- "#555555"

box <- function(x, y, w, h, label, fill = "white", border = mid,
                fontsize = 8.5, fontface = "plain", lineheight = 0.95) {
  grid.rect(x = x, y = y, width = w, height = h,
            gp = gpar(fill = fill, col = border, lwd = 1.0))
  grid.text(label, x = x, y = y,
            gp = gpar(col = dark, fontsize = fontsize, fontface = fontface,
                      lineheight = lineheight))
}

arrow_down <- function(x, y1, y2) {
  grid.lines(x = unit(c(x, x), "npc"), y = unit(c(y1, y2), "npc"),
             arrow = arrow(type = "closed", length = unit(2.2, "mm")),
             gp = gpar(col = mid, lwd = 1.2))
}

draw_figure <- function() {
  grid.newpage()

  # Panel labels
  grid.text("a", 0.035, 0.965, just = c("left", "top"),
            gp = gpar(fontsize = 18, fontface = "bold", col = dark))
  grid.text("b", 0.445, 0.965, just = c("left", "top"),
            gp = gpar(fontsize = 18, fontface = "bold", col = dark))
  grid.text("c", 0.455, 0.485, just = c("left", "top"),
            gp = gpar(fontsize = 18, fontface = "bold", col = dark))

  # Panel A: evidence workflow
  grid.text("Evidence workflow", 0.225, 0.94,
            gp = gpar(fontsize = 10.5, fontface = "bold", col = dark))

  ys <- c(0.865, 0.735, 0.605, 0.475, 0.345, 0.215)
  box(0.225, ys[1], 0.35, 0.075,
      "Genetic resources\nCDD GWAS | hypertension GWAS | cis-eQTL",
      fill = grey_fill, fontsize = 8.7)
  box(0.225, ys[2], 0.29, 0.070,
      "Exploratory shared-SMR screen\nhypothesis-generating",
      fill = "white", fontsize = 8.7)
  box(0.225, ys[3], 0.33, 0.075,
      "Archived four-transcript candidate set\nreconstructed from preserved records",
      fill = blue_fill, border = blue_border, fontsize = 8.3)
  box(0.225, ys[4], 0.30, 0.070,
      "Formal regional characterization\nSMR/HEIDI + pairwise colocalization",
      fill = "white", fontsize = 8.3)
  box(0.225, ys[5], 0.35, 0.075,
      "External and transcriptomic context\nMVP SBP | UKB M50 | bulk | scRNA-seq",
      fill = grey_fill, fontsize = 8.3)
  box(0.225, ys[6], 0.28, 0.065,
      "Qualitative integration\nwithout additive scoring",
      fill = "white", fontsize = 8.5)
  for (i in seq_len(length(ys) - 1)) arrow_down(0.225, ys[i] - 0.047, ys[i + 1] + 0.047)

  grid.text(
    "Historical records do not independently establish the prospective\ntiming of every threshold or analytical component.",
    0.225, 0.105,
    gp = gpar(fontsize = 7.3, fontface = "italic", col = mid, lineheight = 0.95)
  )

  # Panel B: two loci and transcript hypotheses
  grid.text("Candidate loci and transcript hypotheses", 0.745, 0.94,
            gp = gpar(fontsize = 10.5, fontface = "bold", col = dark))

  grid.rect(0.725, 0.78, width = 0.48, height = 0.235,
            gp = gpar(fill = "#EEF5F9", col = blue_border, lwd = 1.5))
  grid.text("Chromosome 10 locus", 0.50, 0.875, just = "left",
            gp = gpar(fontsize = 9.5, fontface = "bold", col = dark))
  box(0.725, 0.835, 0.36, 0.052,
      "NT5C2 | locus-specific follow-up hypothesis",
      fill = blue_fill, border = mid, fontsize = 8.3)
  box(0.725, 0.770, 0.36, 0.052,
      "CNNM2 | alternative transcript hypothesis",
      fill = "white", border = mid, fontsize = 8.3)
  box(0.725, 0.705, 0.36, 0.052,
      "MARCKSL1P1 | alternative transcript hypothesis",
      fill = "white", border = mid, fontsize = 8.3)

  grid.rect(0.725, 0.585, width = 0.48, height = 0.105,
            gp = gpar(fill = brown_fill, col = brown_border, lwd = 1.5))
  grid.text("Chromosome 18 locus", 0.50, 0.62, just = "left",
            gp = gpar(fontsize = 9.5, fontface = "bold", col = dark))
  box(0.725, 0.565, 0.36, 0.052,
      "RMC1 | locus-specific follow-up hypothesis",
      fill = blue_fill, border = mid, fontsize = 8.3)

  grid.text("4 transcript-level hypotheses from 2 independent loci",
            0.725, 0.505,
            gp = gpar(fontsize = 8.8, fontface = "bold", col = dark))

  # Panel C: separate inferential levels
  grid.text("Inference levels", 0.735, 0.465,
            gp = gpar(fontsize = 10.5, fontface = "bold", col = dark))
  grid.lines(x = unit(c(0.72, 0.72), "npc"), y = unit(c(0.135, 0.415), "npc"),
             gp = gpar(col = "#B0B0B0", lwd = 1.2))

  grid.text("Locus-level\ncardiovascular support", 0.585, 0.400,
            gp = gpar(fontsize = 8.7, fontface = "bold", col = dark, lineheight = 0.95))
  grid.text("Transcript-level hypotheses\nremain unresolved", 0.855, 0.400,
            gp = gpar(fontsize = 8.7, fontface = "bold", col = dark, lineheight = 0.95))

  box(0.585, 0.310, 0.22, 0.090,
      "Chromosome 10 locus\nMVP SBP regional support\nnot transcript-specific",
      fill = grey_fill, fontsize = 8.2)
  box(0.585, 0.190, 0.22, 0.090,
      "Chromosome 18 locus\nMVP SBP regional support\nrelated phenotype context",
      fill = grey_fill, fontsize = 8.2)

  box(0.855, 0.330, 0.235, 0.050,
      "NT5C2 | follow-up hypothesis", fill = blue_fill, fontsize = 7.8)
  box(0.855, 0.265, 0.235, 0.050,
      "CNNM2 | alternative hypothesis", fill = "white", fontsize = 7.8)
  box(0.855, 0.200, 0.235, 0.060,
      "MARCKSL1P1 |\nalternative hypothesis", fill = "white", fontsize = 7.0)
  box(0.855, 0.135, 0.235, 0.050,
      "RMC1 | follow-up hypothesis", fill = blue_fill, fontsize = 7.8)

  grid.text("Locus support and transcript attribution were interpreted separately.",
            0.735, 0.070,
            gp = gpar(fontsize = 8.0, fontface = "bold", col = dark))
}

pdf_file <- file.path(out_dir, "Figure1_study_design_FINAL.pdf")
tiff_file <- file.path(out_dir, "Figure1_study_design_FINAL_600dpi.tiff")
png_file <- file.path(out_dir, "Figure1_study_design_FINAL_preview.png")

cairo_pdf(pdf_file, width = 180 / 25.4, height = 140 / 25.4,
          family = "Arial", fallback_resolution = 600)
draw_figure()
dev.off()

ragg::agg_tiff(tiff_file, width = 180, height = 140, units = "mm",
               res = 600, compression = "lzw", background = "white")
draw_figure()
dev.off()

ragg::agg_png(png_file, width = 180, height = 140, units = "mm",
              res = 180, background = "white")
draw_figure()
dev.off()

message("Created:")
message(pdf_file)
message(tiff_file)
message(png_file)
