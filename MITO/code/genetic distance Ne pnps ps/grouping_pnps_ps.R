# pnps_vs_ps_species.R
# Scatter of mean pN/pS vs mean synonymous diversity (pS) with proper Paramecium names + your colors

library(ggplot2)

# -----------------------------
# Data (as you provided)
# -----------------------------
df <- data.frame(
  Species = c(
    "biaurelia","bursaria","caudatum","fokini","multi","pquin",
    "primaurelia","sexaurelia","tetraurelia","triaurelia","tribursaria"
  ),
  Mean_pN_pS = c(
    0.3593980, 0.1330528, 0.1048399, 0.2375698, 0.0762123, 0.4423753,
    0.0401629, 0.0636540, 0.1144409, 0.0895303, 0.1974634
  ),
  Mean_pS = c(
    0.0105606, 0.0275190, 0.1799274, 0.0709674, 0.0116931, 0.0042360,
    0.1680497, 0.2356418, 0.0205187, 0.0274407, 0.0200260
  )
)

# Spearman rho (as you reported)
rho <- -0.6181818

# -----------------------------
# Proper species labels (PNAS-friendly)
# -----------------------------
species_labels <- c(
  biaurelia    = "P. biaurelia",
  bursaria     = "P. bursaria",
  caudatum     = "P. caudatum",
  fokini       = "P. fokini",
  multi        = "P. multimicronucleatum",
  pquin        = "P. quindecaurelia",
  primaurelia  = "P. primaurelia",
  sexaurelia   = "P. sexaurelia",
  tetraurelia  = "P. tetraurelia",
  triaurelia   = "P. triaurelia",
  tribursaria  = "P. tribursaria"
)

# Ensure consistent legend order (optional but nice)
df$Species <- factor(df$Species, levels = names(species_labels))

# -----------------------------
# Your custom color map (mapped to shorthand keys)
# -----------------------------
color_map <- c(
  tribursaria  = "#88c96b",
  bursaria     = "#cc0000",
  multi        = "#9874f4",
  fokini       = "#ffe599",
  caudatum     = "#e69138",
  sexaurelia   = "#559738",
  tetraurelia  = "#0b5394",
  tredecaurelia= "#f1c232",  # not used in this df, kept for consistency
  quadecaurelia= "#ccd9c7",  # not used
  Psch         = "#aaaaaa",  # not used
  pquin        = "#d5ed8a",
  biaurelia    = "#ee79ac",
  triaurelia   = "#d5a6bd",
  primaurelia  = "#6fa8dc"
)

# Subset colors to only those species present in df
color_map_use <- color_map[names(color_map) %in% levels(df$Species)]

# -----------------------------
# Plot
# -----------------------------
p <- ggplot(df, aes(x = Mean_pS, y = Mean_pN_pS, color = Species)) +
  geom_point(size = 3) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linetype = "dashed",
    color = "black",
    linewidth = 0.6
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = expression("Mean synonymous diversity (" * italic(p)[S] * ")"),
    y = expression("Mean " * italic(p)[N] * "/" * italic(p)[S]),
    title = "Negative association between pN/pS and synonymous diversity",
    color = "Species"
  ) +
  scale_color_manual(
    values = color_map_use,
    labels = species_labels
  ) +
  annotate(
    "text",
    x = max(df$Mean_pS) * 0.55,
    y = max(df$Mean_pN_pS) * 0.95,
    label = paste0("Spearman \u03C1 = ", round(rho, 3)),
    size = 5
  )

print(p)

# -----------------------------
# Save
# -----------------------------
ggsave("pnps_vs_ps_species.pdf", p, width = 5.5, height = 4.5)
ggsave("pnps_vs_ps_species.png", p, width = 5.5, height = 4.5, dpi = 300)
