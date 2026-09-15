#!/usr/bin/env Rscript

# Date: 9/9/2026
# Author: Jaaved Mohammed

# Step 08: fig_final_comparison.pdf / .png
#
# Reads results/final_comparison_values.tsv (eight bars) and
# results/final_comparison_nested.tsv (four brackets). Two panels: phastCons470way
# and LINSIGHT. In each panel two x groups (distal peaks, motif occurrences), each
# with four bars: parent set, its background, K27ac subset, its background.
#
# The PDF uses the base pdf() device so every label remains an editable text string.

suppressMessages({
    library(data.table)
    library(ggplot2)
    library(patchwork)
})

RES <- Sys.getenv("RES")
FIG <- Sys.getenv("FIG")

V    <- fread(file.path(RES, "final_comparison_values.tsv"))
NEST <- fread(file.path(RES, "final_comparison_nested.tsv"))


# =============================================================================
# 1. Labels, colours and bar positions
# =============================================================================

GRP <- c("peaks", "k27_peaks", "motifs", "k27_motifs")

LONG <- c(peaks      = "All SNAI2 distal peaks",
          k27_peaks  = "Peaks at H3K27ac-gain enhancers",
          motifs     = "SNAI2 motif occurrences, FIMO p < 1e-4",
          k27_motifs = "Motif occurrences at H3K27ac-gain enhancers")

# legend entries carry the element count, e.g. "All SNAI2 distal peaks (N = 3,877)"
nfg <- V[panel == "phastCons470way"][match(GRP, group), n_fg]
LEG <- setNames(sprintf("%s (N = %s)", LONG[GRP], trimws(format(nfg, big.mark = ","))), GRP)

BLU <- "Matched background"

COL <- c(setNames(c("#b2182b", "#2e8b57", "#6a3d9a", "#e6ab02"), LEG[GRP]),
         setNames("#6f9bc4", BLU))

# x position of each foreground bar; its background bar sits 0.20 to the right
XPOS  <- c(peaks = 0.70, k27_peaks = 1.10, motifs = 1.70, k27_motifs = 2.10)
BAR_W <- 0.185
GREY  <- "grey35"


# =============================================================================
# 2. Helpers
# =============================================================================

# p-value label for a bracket
pfmt <- function(p) {
    ifelse(p == 0,      "p < 1e-300",
    ifelse(p < 0.001,   paste0("p = ", sub("e-0", "e-", formatC(p, format = "e", digits = 0))),
                        paste0("p = ", formatC(p, format = "f", digits = 3))))
}

th <- theme_classic(base_size = 9) +
    theme(axis.text       = element_text(colour = "grey20", size = 8),
          axis.title      = element_text(size = 9),
          axis.ticks.x    = element_blank(),
          plot.title      = element_text(size = 9.5, face = "bold", hjust = 0),
          plot.subtitle   = element_text(size = 8, colour = GREY),
          plot.tag        = element_text(size = 11, face = "bold"),
          legend.title    = element_blank(),
          legend.text     = element_text(size = 7.6),
          legend.key.size = unit(9, "pt"),
          plot.margin     = margin(6, 12, 4, 6))


# =============================================================================
# 3. One panel
# =============================================================================

panel_plot <- function(panel_, ylab, title, sub) {

    # -- the eight bars of this panel ----------------------------------------
    v <- V[panel == panel_][match(GRP, group)]
    v[, x  := XPOS[group]]           # foreground bar
    v[, xb := XPOS[group] + 0.20]    # its background bar

    fg <- v[, .(xpos = x,  who = LEG[group], value = obs,     sd = NA_real_)]
    bg <- v[, .(xpos = xb, who = BLU,        value = bg_mean, sd = bg_sd)]

    d <- rbind(fg, bg)
    d[, who := factor(who, levels = c(LEG[GRP], BLU))]

    ymax <- max(d$value + fifelse(is.na(d$sd), 0, d$sd))

    # -- lower brackets: each foreground bar vs its own background ---------------
    br <- v[, .(x1  = x,
                x2  = xb,
                lab = pfmt(p),
                y   = pmax(obs, bg_mean + bg_sd) + ymax * 0.05)]

    # -- upper brackets: subset vs parent set (nested test) ----------------------
    nb <- data.table(x1  = XPOS[c("peaks", "motifs")],
                     x2  = XPOS[c("k27_peaks", "k27_motifs")],
                     lab = pfmt(NEST[panel == panel_, p]),
                     y   = c(max(br$y[1:2]), max(br$y[3:4])) + ymax * 0.10)

    bars <- rbind(br, nb)

    # short vertical ticks at both ends of every bracket
    tick <- rbind(bars[, .(x = x1, y = y - ymax * 0.016, ye = y)],
                  bars[, .(x = x2, y = y - ymax * 0.016, ye = y)])

    # -- draw ----------------------------------------------------------------------
    ggplot(d, aes(xpos, value, fill = who)) +
        geom_col(width = BAR_W) +
        geom_errorbar(aes(ymin = value - sd, ymax = value + sd),
                      width = 0.07, linewidth = 0.35, colour = GREY, na.rm = TRUE) +
        geom_segment(data = bars, aes(x = x1, xend = x2, y = y, yend = y),
                     inherit.aes = FALSE, linewidth = 0.3, colour = GREY) +
        geom_segment(data = tick, aes(x = x, xend = x, y = y, yend = ye),
                     inherit.aes = FALSE, linewidth = 0.3, colour = GREY) +
        geom_text(data = bars, aes(x = (x1 + x2) / 2, y = y, label = lab),
                  inherit.aes = FALSE, vjust = -0.45, size = 2.4, colour = GREY) +
        scale_fill_manual(values = COL, drop = FALSE,
                          guide = guide_legend(nrow = 2, byrow = TRUE)) +
        scale_x_continuous(breaks = 1:2, limits = c(0.55, 2.45),
                           labels = c("Distal peaks", "Motif occurrences")) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.26))) +
        labs(x = NULL, y = ylab, title = title, subtitle = sub) +
        th
}


# =============================================================================
# 4. Assemble and write
# =============================================================================

pa <- panel_plot("phastCons470way",
                 "Fraction of aligned bases with phastCons470way >= 0.8",
                 "Cross-species conservation", "470-mammal alignment")

pb <- panel_plot("LINSIGHT",
                 "Mean LINSIGHT score",
                 "Human polymorphism and divergence", "Probability of negative selection")

fig <- (pa | pb) +
    plot_annotation(tag_levels = "a") +
    plot_layout(guides = "collect") &
    theme(legend.position = "top")

pdf(file.path(FIG, "fig_final_comparison.pdf"), width = 10.5, height = 5.0, useDingbats = FALSE)
print(fig)
invisible(dev.off())

ggsave(file.path(FIG, "fig_final_comparison.png"), fig, width = 10.5, height = 5.0, dpi = 300, bg = "white")

cat("wrote", file.path(FIG, "fig_final_comparison.pdf"), "\n")
