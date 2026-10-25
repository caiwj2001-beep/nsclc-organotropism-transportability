# =============================================================================
# Figures 1-3 | NSCLC organotropism manuscript V1
# 规则：source(fig_helpers.R)；框宽由 box_size_from_text 计算；流程排位用 flow_layout；
#       森林图注释列用 forest_annot_table；每张图跑 detect_overlap；白底、300 dpi。
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(grid); library(gridExtra); library(dplyr); library(scales)
})
source("E:/科研文件/fig_helpers.R")

PROJ <- "E:/科研文件/肺癌器官特异性转移"
FD   <- file.path(PROJ, "05_results/figdata")
FIG  <- file.path(PROJ, "04_figures"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

NAVY <- "#1F3A5F"; CORAL <- "#E4776E"; GREY <- "#8A8F98"
HCOL <- c(LUAD = NAVY, LUSC = CORAL)
theme_set(theme_classic(base_size = 9) +
            theme(plot.background = element_rect(fill = "white", colour = NA),
                  panel.background = element_rect(fill = "white", colour = NA),
                  strip.background = element_rect(fill = "white", colour = NA),
                  axis.line = element_line(colour = "black", linewidth = 0.35),
                  axis.ticks = element_line(colour = "black", linewidth = 0.35),
                  legend.key = element_blank(), plot.margin = margin(8, 8, 8, 8)))

rd <- function(f) read.csv(file.path(FD, f), check.names = FALSE)
ovl <- list()

# ---------------------------------------------------------------- Figure 1 ---
# A. flow
flow_labels <- c(
  "MSK-CHORD NSCLC samples sequenced\nn = 7,809",
  "Patients mapped to a mutation record\nn = 7,466",
  "Main analysis set\nn = 6,516  (LUAD 5,717 | LUSC 799)")
fl <- flow_layout(flow_labels, direction = "down", spacing_in = 0.42,
                  fontsize = 8, plot_w_in = 3.5, plot_h_in = 3.4, x_center_norm = 40)

# arithmetic guard: every excluded patient must be accounted for
stopifnot(7809 - 343 == 7466)
stopifnot(7466 - 250 - 148 - 3 - 473 - 76 == 6516)

exc <- c("No mutation record
n = 343",
         "Neuroendocrine
n = 250",
         "Other histology
n = 148",
         "Unclassifiable
n = 3",
         "NSCLC-NOS
(descriptive only) n = 473",
         "Poorly differentiated
n = 76")
ex_sz <- lapply(exc, box_size_from_text, fontsize = 6.0,
                text_pad_in = 0.055, line_height_in = 0.125)
ex_h  <- sapply(ex_sz, function(s) s["height_in"])
ex_gap <- 0.045
ex_top <- 93
ex_y <- numeric(length(exc)); run <- 0
for (i in seq_along(exc)) { ex_y[i] <- ex_top - to_norm(run + ex_h[i] / 2, 3.4); run <- run + ex_h[i] + ex_gap }
ex_w <- to_norm(sapply(ex_sz, function(s) s["width_in"]), 3.5)

pA <- ggplot() +
  geom_rect(data = fl, aes(xmin = x - w/2, xmax = x + w/2, ymin = y - h/2, ymax = y + h/2),
            fill = "white", colour = NAVY, linewidth = 0.45) +
  geom_text(data = fl, aes(x = x, y = y, label = label), size = 2.55, lineheight = 0.95) +
  geom_segment(data = data.frame(x = fl$x[1:2], y = fl$y[1:2] - fl$h[1:2]/2,
                                 xend = fl$x[2:3], yend = fl$y[2:3] + fl$h[2:3]/2),
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.06, "in")), linewidth = 0.35) +
  geom_rect(data = data.frame(x = 84, y = ex_y, w = ex_w, h = to_norm(ex_h, 3.4)),
            aes(xmin = x - w/2, xmax = x + w/2, ymin = y - h/2, ymax = y + h/2),
            fill = "white", colour = GREY, linewidth = 0.3, linetype = "22") +
  geom_text(data = data.frame(x = 84, y = ex_y, label = exc), aes(x = x, y = y, label = label),
            size = 1.95, colour = "#4A4F57", lineheight = 0.88) +
  geom_segment(aes(x = fl$x[1] + fl$w[1]/2, xend = 84 - ex_w[1]/2,
                   y = fl$y[1], yend = ex_y[1]), linewidth = 0.28, linetype = "22", colour = GREY) +
  geom_segment(aes(x = fl$x[2] + fl$w[2]/2, xend = 84 - ex_w[2]/2,
                   y = fl$y[2], yend = (ex_y[2] + ex_y[6])/2), linewidth = 0.28, linetype = "22", colour = GREY) +
  coord_cartesian(xlim = c(0, 100), ylim = c(2, 100), expand = FALSE) +
  theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA),
                       plot.margin = margin(12, 12, 12, 18)) +
  labs(tag = "A") + theme(plot.tag = element_text(size = 11, face = "bold"))

# B. organ distribution among distant-metastatic patients
dB <- rd("fig1b_organ.csv") %>% tidyr::pivot_longer(-organ, names_to = "hist", values_to = "pct")
pB <- ggplot(dB, aes(x = factor(organ, levels = c("Bone","CNS","Pleura","Liver","Adrenal")),
                     y = pct, fill = hist)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.66, colour = NA) +
  geom_text(aes(label = sprintf("%.1f", pct)), position = position_dodge(width = 0.72),
            vjust = -0.35, size = 2.1, show.legend = FALSE) +
  scale_fill_manual(values = HCOL, name = NULL) +
  scale_y_continuous(limits = c(0, 78), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Patients with organ involvement (%)", tag = "B") +
  theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
        legend.background = element_rect(fill = alpha("white", 0.75), colour = NA),
        legend.key.size = unit(0.28, "cm"), legend.text = element_text(size = 7),
        plot.tag = element_text(size = 11, face = "bold"))

# C. single-organ distribution
dC <- rd("fig1c_single_organ.csv") %>% tidyr::pivot_longer(-organ, names_to = "hist", values_to = "pct")
pC <- ggplot(dC, aes(x = factor(organ, levels = c("Bone","CNS","Pleura","Liver","Adrenal")),
                     y = pct, fill = hist)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.66, colour = NA) +
  geom_text(aes(label = sprintf("%.1f", pct)), position = position_dodge(width = 0.72),
            vjust = -0.35, size = 2.1, show.legend = FALSE) +
  scale_fill_manual(values = HCOL, name = NULL) +
  scale_y_continuous(limits = c(0, 30), expand = expansion(mult = c(0, 0.06))) +
  labs(x = NULL, y = "Single-organ patients (%)", tag = "C",
       subtitle = expression(chi^2 * " = 4.70, df = 6, " * italic(P) * " = .589")) +
  theme(legend.position = "none", plot.subtitle = element_text(size = 7, colour = "#4A4F57"),
        plot.tag = element_text(size = 11, face = "bold"))

# D. mutation landscape (dumbbell: two rates per gene, no end labels -> no collisions)
dD <- rd("fig1d_mutation.csv")
ordD <- dD %>% arrange(pmax(LUAD, LUSC)) %>% pull(gene)
dD$gene <- factor(dD$gene, levels = ordD)
pD <- ggplot(dD, aes(y = gene)) +
  geom_segment(aes(x = LUAD, xend = LUSC, yend = gene), colour = "#C9CED6", linewidth = 0.5) +
  geom_point(aes(x = LUAD, fill = "LUAD"), size = 2.15, shape = 21, colour = "white", stroke = 0.35) +
  geom_point(aes(x = LUSC, fill = "LUSC"), size = 1.45, shape = 22, colour = "white", stroke = 0.35) +
  scale_fill_manual(values = HCOL, name = NULL) +
  scale_x_continuous(limits = c(0, 95), breaks = c(0, 25, 50, 75),
                     expand = expansion(mult = c(0.01, 0.02))) +
  labs(x = "Mutated (%)", y = NULL, tag = "D") +
  theme(legend.position = c(0.97, 0.06), legend.justification = c(1, 0),
        legend.background = element_rect(fill = alpha("white", 0.8), colour = NA),
        legend.key.size = unit(0.28, "cm"), legend.text = element_text(size = 6.6),
        axis.text.y = element_text(size = 6.6, face = "italic"),
        plot.tag = element_text(size = 11, face = "bold"))

# E. CNA landscape (dumbbell, faceted by alteration type)
dE <- rd("fig1e_cna.csv")
ordE <- dE %>% arrange(type, pmax(LUAD, LUSC)) %>% pull(gene)
dE$gene <- factor(dE$gene, levels = ordE)
pE <- ggplot(dE, aes(y = gene)) +
  geom_segment(aes(x = LUAD, xend = LUSC, yend = gene), colour = "#C9CED6", linewidth = 0.5) +
  geom_point(aes(x = LUAD, fill = "LUAD"), size = 2.0, shape = 21, colour = "white", stroke = 0.35) +
  geom_point(aes(x = LUSC, fill = "LUSC"), size = 1.4, shape = 22, colour = "white", stroke = 0.35) +
  scale_fill_manual(values = HCOL, name = NULL) +
  facet_grid(type ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(limits = c(0, 21), breaks = c(0, 5, 10, 15),
                     expand = expansion(mult = c(0.02, 0.03))) +
  labs(x = "Altered (%)", y = NULL, tag = "E") +
  theme(legend.position = "none", strip.text = element_text(size = 6.6),
        axis.text.y = element_text(size = 6.2, face = "italic"),
        plot.tag = element_text(size = 11, face = "bold"))

for (nm in c("A","B","C","D","E")) ovl[[paste0("Fig1", nm)]] <- detect_overlap(get(paste0("p", nm)))

fig1 <- arrangeGrob(pA, pB, pC, pD, pE, ncol = 2,
                    layout_matrix = rbind(c(1, 2), c(3, 4), c(5, 5)),
                    heights = c(1.05, 0.85, 0.95))
ggsave(file.path(FIG, "Figure1.pdf"), fig1, width = 7.2, height = 8.6, device = cairo_pdf, bg = "white")
ggsave(file.path(FIG, "Figure1.tiff"), fig1, width = 7.2, height = 8.6, dpi = 300, compression = "lzw", bg = "white")
ggsave(file.path(FIG, "Figure1.png"), fig1, width = 7.2, height = 8.6, dpi = 300, bg = "white")
cat("Figure1 done\n")

# ---------------------------------------------------------------- Figure 2 ---
f2 <- rd("fig2_forest.csv")
lvl <- unique(f2$label)                      # file order = Figure 3 order (adrenal first)
f2$label <- factor(f2$label, levels = rev(lvl))
f2$y <- as.numeric(f2$label)
lu <- subset(f2, stratum == "LUAD"); ls_ <- subset(f2, stratum == "LUSC")
# match by label (never by position)
grab <- function(d, col) d[[col]][match(levels(f2$label), as.character(d$label))]
lu <- data.frame(label = levels(f2$label), y = as.numeric(levels(f2$label)) * 0 + seq_along(levels(f2$label)),
                 OR = grab(lu, "OR"), lo = grab(lu, "lo"), hi = grab(lu, "hi"),
                 p = grab(lu, "p"), p_int = grab(lu, "p_int"), mde = grab(lu, "mde"))
ls_ <- data.frame(label = levels(f2$label),
                  OR = grab(ls_, "OR"), lo = grab(ls_, "lo"), hi = grab(ls_, "hi"))
lu$y <- (length(levels(f2$label)):1)
ls_$y <- lu$y
YL <- c(0.2, length(lvl) + 1.6); HEAD <- length(lvl) + 0.9

gg_f <- ggplot() +
  geom_vline(xintercept = 1, linewidth = 0.3, colour = GREY, linetype = "22") +
  geom_errorbar(data = lu, aes(y = y + 0.17, xmin = lo, xmax = hi), orientation = "y",
                width = 0.11, linewidth = 0.5, colour = NAVY) +
  geom_errorbar(data = ls_, aes(y = y - 0.17, xmin = lo, xmax = hi), orientation = "y",
                width = 0.11, linewidth = 0.5, colour = CORAL) +
  geom_point(data = lu, aes(y = y + 0.17, x = OR), shape = 21, size = 1.9, fill = "white",
             colour = NAVY, stroke = 0.7) +
  geom_point(data = ls_, aes(y = y - 0.17, x = OR), shape = 22, size = 1.7, fill = CORAL, colour = CORAL) +
  scale_x_log10(breaks = c(0.5, 1, 2, 5, 10), labels = c("0.5","1","2","5","10"),
                limits = c(0.22, 14)) +
  scale_y_continuous(breaks = lu$y, labels = levels(f2$label), limits = YL, expand = c(0, 0)) +
  labs(x = "Odds ratio (log scale)", y = NULL, tag = "A") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 7), plot.tag = element_text(size = 11, face = "bold"))

# annotation panel: same y scale; contrast label carried in the table (self-identifying rows)
tab <- data.frame(y = c(HEAD, lu$y),
                  c0 = c("Contrast", as.character(lu$label)),
                  c1 = c("LUAD OR (95% CI)", sprintf("%.2f (%.2f\u2013%.2f)", lu$OR, lu$lo, lu$hi)),
                  c2 = c("LUSC OR (95% CI)", sprintf("%.2f (%.2f\u2013%.2f)", ls_$OR, ls_$lo, ls_$hi)),
                  c3 = c("P int.", ifelse(lu$p_int < 0.001, "<.001", sprintf("%.2f", lu$p_int))),
                  c4 = c("MDE", sprintf("%.1f", lu$mde)),
                  stringsAsFactors = FALSE)
XS <- c(0.00, 0.315, 0.575, 0.765, 0.875)
lay <- do.call(rbind, lapply(seq_along(XS), function(k) {
  data.frame(y = tab$y, x = XS[k], lab = tab[[k + 1]], bold = tab$y == HEAD,
             stringsAsFactors = FALSE)
}))
gg_t <- ggplot(lay, aes(x = x, y = y, label = lab)) +
  geom_text(aes(fontface = ifelse(bold, "bold", "plain")), hjust = 0, size = 2.28) +
  scale_y_continuous(limits = YL, expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0.2, 0.5, 0.8), labels = c("0.5", "1", "2")) +
  labs(x = "Odds ratio (log scale)", y = NULL) +
  theme_classic(base_size = 7) +
  theme(axis.text.x = element_text(colour = "white", size = 7),
        axis.ticks.x = element_line(colour = "white"),
        axis.title.x = element_text(colour = "white", size = 7),
        axis.line.x = element_line(colour = "white"),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = margin(t = 8, r = 2, b = 8, l = 2))

lg <- tableGrob(data.frame(x = c("\u25cb LUAD", "\u25a0 LUSC"), check.names = FALSE), rows = NULL,
                theme = ttheme_minimal(base_size = 7,
                                       core = list(fg_params = list(hjust = 0, x = 0.02, col = c(NAVY, CORAL)),
                                                   bg_params = list(fill = NA)),
                                       colhead = list(fg_params = list(col = "white")),
                                       padding = unit(c(0.8, 2), "mm")))

fig2 <- arrangeGrob(gg_f, arrangeGrob(gg_t, lg, ncol = 1, heights = c(5.4, 0.55)),
                    ncol = 2, widths = c(2.05, 4.05))
ggsave(file.path(FIG, "Figure2.pdf"), fig2, width = 7.4, height = 3.9, device = cairo_pdf, bg = "white")
ggsave(file.path(FIG, "Figure2.tiff"), fig2, width = 7.4, height = 3.9, dpi = 300, compression = "lzw", bg = "white")
ggsave(file.path(FIG, "Figure2.png"), fig2, width = 7.4, height = 3.9, dpi = 300, bg = "white")
cat("Figure2 done\n")

# ---------------------------------------------------------------- Figure 3 ---
d3a <- rd("fig3a_scenarios.csv")
d3a$scenario_en <- factor(d3a$scenario_en, levels = c("Complete case","Missing indicator","Primary sample","Single organ"))
p3a <- ggplot(d3a, aes(x = OR, y = scenario_en, colour = label)) +
  geom_vline(xintercept = 1, linewidth = 0.3, colour = GREY, linetype = "22") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.14, linewidth = 0.5) +
  geom_point(size = 1.9) +
  facet_wrap(~ label, ncol = 1, strip.position = "top") +
  scale_x_log10(breaks = c(0.4, 0.6, 1, 1.5, 2.5, 4), limits = c(0.3, 4.2)) +
  scale_colour_manual(values = c("STK11 \u2013 adrenal" = NAVY, "SMARCA4 \u2013 adrenal" = CORAL), guide = "none") +
  labs(x = "Adjusted odds ratio", y = NULL, tag = "A") +
  theme(strip.text = element_text(size = 7), axis.text.y = element_text(size = 6.6),
        axis.text.x = element_text(size = 6.8), plot.tag = element_text(size = 11, face = "bold"))

d3b <- rd("fig3b_cooccurrence.csv")
d3b$lab <- sprintf("%s\nOR %.1f, n = %d co-altered", d3b$pair, d3b$OR, d3b$both)
p3b <- ggplot(d3b, aes(x = OR, y = reorder(pair, OR))) +
  geom_col(width = 0.55, fill = NAVY, alpha = 0.85) +
  geom_text(aes(label = sprintf("OR %.1f (n = %d)", OR, both)), hjust = -0.12, size = 2.3) +
  scale_x_continuous(limits = c(0, 17), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Co-occurrence odds ratio", y = NULL, tag = "B") +
  theme(axis.text.y = element_text(size = 6.8), plot.tag = element_text(size = 11, face = "bold"))

d3c <- rd("fig3c_attenuation.csv")
d3c$model <- factor(d3c$model, levels = d3c$model)
p3c <- ggplot(d3c, aes(x = model, y = OR)) +
  geom_hline(yintercept = 1, linewidth = 0.3, colour = GREY, linetype = "22") +
  geom_line(aes(group = 1), linewidth = 0.45, colour = GREY) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.09, linewidth = 0.45, colour = NAVY) +
  geom_point(size = 2, shape = 21, fill = "white", colour = NAVY, stroke = 0.75) +
  geom_text(aes(label = sprintf("%.2f", OR)), vjust = -1.15, size = 2.2) +
  scale_y_continuous(limits = c(0.85, 2.15), breaks = c(1, 1.5, 2)) +
  labs(x = NULL, y = "SMARCA4 \u2013 bone OR", tag = "C") +
  theme(axis.text.x = element_text(size = 6.2, angle = 18, hjust = 0.9),
        axis.text.y = element_text(size = 6.8), plot.tag = element_text(size = 11, face = "bold"))

fig3 <- arrangeGrob(p3a, p3b, p3c, ncol = 2, layout_matrix = rbind(c(1, 1), c(2, 3)),
                    heights = c(1.25, 1))
ggsave(file.path(FIG, "Figure3.pdf"), fig3, width = 7.2, height = 6.0, device = cairo_pdf, bg = "white")
ggsave(file.path(FIG, "Figure3.tiff"), fig3, width = 7.2, height = 6.0, dpi = 300, compression = "lzw", bg = "white")
ggsave(file.path(FIG, "Figure3.png"), fig3, width = 7.2, height = 6.0, dpi = 300, bg = "white")
ovl[["Fig3A"]] <- detect_overlap(p3a); ovl[["Fig3C"]] <- detect_overlap(p3c)
cat("Figure3 done\n")

cat("\n=== overlap check (NULL = clean) ===\n")
for (nm in names(ovl)) {
  cat(nm, ": ")
  if (is.null(ovl[[nm]]) || (is.data.frame(ovl[[nm]]) && nrow(ovl[[nm]]) == 0)) cat("clean\n") else {
    cat(nrow(ovl[[nm]]), "overlap(s)\n"); print(utils::head(ovl[[nm]], 5))
  }
}
