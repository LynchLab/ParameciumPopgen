indr  <- 'genetic_vs_geographic/'
fname <- '*_genetic_vs_geographic.csv'
fpl   <- 'genetic_vs_geographic/species_geo_gen_cor.pdf'

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tools)
})

# ---- Read & combine all species files ----
files <- list.files(indr, full.names = FALSE)
stopifnot(length(files) > 0)

df_list <- lapply(files, function(f) {
  d <- read.csv(file.path(indr, f))
  d$species <- file_path_sans_ext(basename(f))
  d
})

df_all <- bind_rows(df_list)

# ---- Print unique species names ----
cat("Unique species names found:\n")
print(unique(df_all$species))

# ---- Per-species scaling & filter ----
df_scaled <- df_all %>%
  group_by(species) %>%
  mutate(Genetic_Distance = Genetic_Distance / max(Genetic_Distance, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(species) %>%
  filter(n() > 1) %>%
  ungroup()

# Keep facet order matching input file order
species_order <- file_path_sans_ext(files)
species_order <- species_order[species_order %in% unique(df_scaled$species)]

# ---- Define color mapping for species ----
species_colors <- c(
  "bursaria_genetic_vs_geographic"      = "#f76666",  # burs
  "multi_genetic_vs_geographic"         = "#9874f4",  # mult
  "fokini_genetic_vs_geographic"        = "#ffe599",  # foki
  "caudatum_genetic_vs_geographic"      = "#e69138",  # caud
  "sexaurelia_genetic_vs_geographic"    = "#559738",  # sexa
  "biaurelia_genetic_vs_geographic"     = "#ee79ac",  # biau
  "primaurelia_genetic_vs_geographic"   = "#6fa8dc",  # prima
  "tetraurelia_genetic_vs_geographic"   = "#0b5394",  # tetr
  "triaurelia_genetic_vs_geographic"    = "#d5a6bd",  # tria
  "tredecaurelia_genetic_vs_geographic" = "#f1c232",  # tred
  "tribursaria_genetic_vs_geographic"   = "#ff7f00",  # tribursaria (custom color)
  "pquin_genetic_vs_geographic"         = "#999999"   # pquin (neutral gray)
)

# ---- Spearman correlation per species ----
cors_tbl <- df_scaled %>%
  group_by(species) %>%
  summarise(
    rho    = suppressWarnings(cor(Geographic_Distance, Genetic_Distance, method = "spearman")),
    pvalue = suppressWarnings(cor.test(Geographic_Distance, Genetic_Distance, method = "spearman")$p.value),
    .groups = "drop"
  )

# Significant positive correlations
cors_sigpos <- cors_tbl %>% filter(rho > 0, pvalue < 0.05)

# Create panel labels
cors_tbl <- cors_tbl %>%
  mutate(strip_lab = sprintf("%s  %.2f%s", species, rho, ifelse(pvalue < 0.05, "**", "")))

# Join back to main dataframe
df_plot <- df_scaled %>%
  left_join(cors_tbl %>% select(species, strip_lab), by = "species") %>%
  mutate(
    strip_lab = factor(strip_lab,
                       levels = cors_tbl$strip_lab[match(species_order, cors_tbl$species)])
  )




# ---- Plot ----
p <- ggplot(df_plot, aes(x = Geographic_Distance, y = Genetic_Distance, color = species)) +
  geom_point(alpha = 0.7, size = 0.8) +
  facet_wrap(~ strip_lab, ncol = 2, nrow = 6, scales = "fixed") +
  scale_x_continuous(limits = c(0, 2e4), breaks = seq(0, 2e4, by = 2e3), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1.1, by = 0.2), expand = c(0, 0)) +
  scale_color_manual(values = species_colors) +
  labs(x = "Geographic distance", y = "Genetic distance", color = "Species") +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90"),
    strip.text = element_text(face = "bold"),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8)),
    legend.position = "bottom",
    legend.text = element_text(size = 8)
  )

# ---- Save PDF ----
ggsave(filename = fpl, plot = p, width = 8.5, height = 11, units = "in")

cat("PDF saved as:", fpl, "\n")





