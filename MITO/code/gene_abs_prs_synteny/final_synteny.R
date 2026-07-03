# ==========================================
# Mitochondrial gene order plot (2 columns)
# Clean version:
# - removes "Left" / "Right"
# - removes top title
# - makes genomic coordinate label bigger
# - keeps bur/trib first ordering internally, but does not print it
# - swaps multi and fokini positions
# ==========================================

library(ggplot2)
library(gggenes)
library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(grid)

# ---------- INPUT / OUTPUT ----------
infile <- "~/UGENE_Data/gene_coordinates.xlsx"
sheet_name <- 1
outpng <- "~/UGENE_Data/gene_order_ALL_species_2cols_clean.png"

# ---------- FIXED GENOME SIZES (bp) ----------
genome_size_bp <- c(
  "biaurelia"          = 39149,
  "decaurelia"         = 40096,
  "dodecaurelia"       = 40002,
  "jenningsi"          = 40071,
  "novaurelia"         = 39603,
  "octaurelia"         = 39795,
  "pentaurelia"        = 39838,
  "primaurelia"        = 39827,
  "psch"               = 39803,
  "sexaurelia"         = 38678,
  "sonneborni"         = 40237,
  "tetraurelia"        = 39984,
  "triaurelia"         = 39599,
  "caudatum"           = 43875,
  "multi"              = 39096,
  "fokini"             = 29946,
  "bursaria"           = 44869,
  "tribursaria"        = 48309
)

# ---------- READ XLSX ----------
gene_data <- read_excel(infile, sheet = sheet_name) %>%
  rename(
    sample_name = species,
    stop = end
  ) %>%
  mutate(
    sample_name = as.character(sample_name),
    gene_name   = as.character(gene_name),
    start       = as.numeric(start),
    stop        = as.numeric(stop),
    strand      = as.character(strand),
    feature     = as.character(feature)
  ) %>%
  filter(!is.na(sample_name), !is.na(gene_name), !is.na(start), !is.na(stop))

# ---------- SPECIES LABELS ----------
species_label_map <- c(
  "primaurelia"        = "P. primaurelia",
  "biaurelia"          = "P. biaurelia",
  "triaurelia"         = "P. triaurelia",
  "tetraurelia"        = "P. tetraurelia",
  "pentaurelia_87"     = "P. pentaurelia",
  "sexaurelia"         = "P. sexaurelia",
  "octaurelia_138"     = "P. octaurelia",
  "octaurelia_139"     = "P. octaurelia",
  "novaurelia_TE"      = "P. novaurelia",
  "decaurelia_223"     = "P. decaurelia",
  "dodecaurelia_274"   = "P. dodecaurelia",
  "tredecaurelia_209"  = "P. tredecaurelia",
  "quadecaurelia"      = "P. quadecaurelia",
  "quindecaurelia"     = "P. quindecaurelia",
  "sonneborni_30995"   = "P. sonneborni",
  "jenningsi_M"        = "P. jenningsi",
  "G4__Psche"          = "P. schweiakoffi",
  "caudatum"           = "P. caudatum",
  "multi"              = "P. multimicronucleatum",
  "pfoki"              = "P. fokini",
  "bursaria"           = "P. bursaria",
  "tribursaria"        = "P. tribursaria"
)

gene_data <- gene_data %>%
  mutate(
    species_label = ifelse(
      sample_name %in% names(species_label_map),
      unname(species_label_map[sample_name]),
      sample_name
    )
  )

# ---------- OPTIONAL: REMOVE UNWANTED ----------
gene_data <- gene_data %>%
  filter(!sample_name %in% c("CL13-6__Pcalkinski", "Pa_a_PFE_42__tredecaurelia"))

# ---------- FLIP RULES ----------
flip_species <- c(
  "G4__Psche", "bursaria", "tribursaria", "quadecaurelia",
  "sexaurelia", "tetraurelia", "triaurelia", "caudatum"
)
flip_samples <- unique(flip_species)

# ---------- ORDER: bursaria & tribursaria first ----------
ordered_samples <- c(
  intersect(c("bursaria", "tribursaria"), unique(gene_data$sample_name)),
  setdiff(unique(gene_data$sample_name), c("bursaria", "tribursaria"))
)

# swap multi and fokini (correct names)
i_multi  <- which(ordered_samples %in% c("multi"))
i_fokini <- which(ordered_samples %in% c("fokini", "pfoki"))

if (length(i_multi) == 1 && length(i_fokini) == 1) {
  tmp <- ordered_samples[i_multi]
  ordered_samples[i_multi] <- ordered_samples[i_fokini]
  ordered_samples[i_fokini] <- tmp
}

# ---------- 2-COLUMN LAYOUT ----------
n_total <- length(ordered_samples)
n_left  <- ceiling(n_total / 2)

panel_df <- tibble(sample_name = ordered_samples) %>%
  mutate(
    panel = factor(ifelse(row_number() <= n_left, "Left", "Right"),
                   levels = c("Left", "Right"))
  ) %>%
  group_by(panel) %>%
  mutate(y = rev(row_number())) %>%
  ungroup()

# ---------- GENOME LENGTHS ----------
genome_lengths <- tibble(sample_name = ordered_samples) %>%
  mutate(genome_length = as.numeric(genome_size_bp[sample_name])) %>%
  left_join(
    gene_data %>%
      filter(sample_name %in% ordered_samples) %>%
      mutate(max_pos = pmax(start, stop, na.rm = TRUE)) %>%
      group_by(sample_name) %>%
      summarise(computed_len = max(max_pos, na.rm = TRUE), .groups = "drop"),
    by = "sample_name"
  ) %>%
  mutate(genome_length = ifelse(is.na(genome_length), computed_len, genome_length)) %>%
  select(sample_name, genome_length)

# ---------- PLOT DATA ----------
gene_plot_data <- gene_data %>%
  filter(sample_name %in% ordered_samples) %>%
  left_join(genome_lengths, by = "sample_name") %>%
  mutate(
    flip = sample_name %in% flip_samples,
    start_raw = pmin(start, stop),
    end_raw   = pmax(start, stop),
    start_adj = ifelse(flip, genome_length - end_raw, start_raw),
    end_adj   = ifelse(flip, genome_length - start_raw, end_raw),
    
    forward = case_when(
      !flip & strand == "+" ~ TRUE,
      !flip & strand == "-" ~ FALSE,
      flip  & strand == "+" ~ FALSE,
      flip  & strand == "-" ~ TRUE,
      TRUE ~ start < stop
    ),
    
    gene_class = case_when(
      str_detect(feature, regex("tRNA", ignore_case = TRUE)) ~ "tRNA",
      str_detect(feature, regex("rRNA", ignore_case = TRUE)) ~ "rRNA",
      str_detect(gene_name, "^R") ~ "Conserved",
      str_detect(gene_name, "^G") ~ "Species_specific",
      TRUE ~ "Protein_coding"
    ),
    gene_class = factor(
      gene_class,
      levels = c("tRNA", "rRNA", "Conserved", "Species_specific", "Protein_coding")
    )
  ) %>%
  left_join(panel_df, by = "sample_name")

# ---------- BASELINES ----------
scale_lines <- genome_lengths %>%
  left_join(panel_df, by = "sample_name") %>%
  mutate(y_scale = y - 0.5)

# ---------- LABELS ----------
label_df <- genome_lengths %>%
  left_join(panel_df, by = "sample_name") %>%
  left_join(gene_data %>% distinct(sample_name, species_label), by = "sample_name") %>%
  mutate(
    label = paste0(species_label, " (", genome_length, " bp)"),
    y_text = y + 0.4,
    x_text = genome_length / 2
  )

# ---------- COLORS ----------
gene_class_colors <- c(
  "tRNA"             = "#d73027",
  "rRNA"             = "#4575b4",
  "Conserved"        = "#1a9850",
  "Species_specific" = "#f781bf",
  "Protein_coding"   = "#fdae61"
)

# ---------- AUTO NG SIZE ----------
png_width_in  <- 13
row_height_in <- 0.9
png_height_in <- max(10, (n_left * row_height_in) + 3)

# ---------- PLOT ----------
p <- ggplot(
  gene_plot_data,
  aes(xmin = start_adj, xmax = end_adj, y = y, fill = gene_class, forward = forward)
) +
  geom_segment(
    data = scale_lines,
    aes(x = 0, xend = genome_length, y = y_scale, yend = y_scale),
    inherit.aes = FALSE,
    linewidth = 0.35,
    color = "white"
  ) +
  geom_gene_arrow(
    arrowhead_height = unit(8, "mm"),
    arrowhead_width  = unit(4, "mm"),
    height           = unit(4, "mm")
  ) +
  facet_wrap(~ panel, ncol = 2, scales = "free_x") +
  geom_text(
    data = label_df,
    aes(x = x_text, y = y_text, label = label),
    inherit.aes = FALSE,
    size = 8
  ) +
  scale_fill_manual(
    values = gene_class_colors,
    breaks = c("tRNA", "rRNA", "Species_specific", "Protein_coding"),
    labels = c(
      "tRNA",
      "rRNA",
      "Species-specific genes",
      "Protein-coding genes"
    ),
    name = "Gene class"
  ) +
  labs(
    title = "Gene arrangements",
    x = "Genomic coordinate (bp)",
    y = NULL
  ) +
  theme_genes() +
  theme(
    strip.text = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.title.x = element_text(size = 24, face = "bold", color = "black"),
    axis.line.x = element_line(color = "black", linewidth = 0.8),
    axis.ticks.x = element_line(color = "black", linewidth = 0.6),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    panel.spacing = unit(2, "lines"),
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(15, 15, 15, 15),
    plot.title = element_blank()
  )

# ---------- SAVE ----------
ggsave(
  filename = outpng,
  plot     = p,
  width    = png_width_in,
  height   = png_height_in,
  units    = "in",
  dpi      = 600,
  bg       = "white"
)

cat("Saved:", outpng, "\n")
cat("Plotted species:", n_total, " (Left:", n_left, ", Right:", n_total - n_left, ")\n")
cat("PNG size (in):", png_width_in, "x", png_height_in, " at 600 dpi\n")

cat("\nGenome sizes used (bp):\n")
print(
  genome_lengths %>%
    left_join(gene_data %>% distinct(sample_name, species_label), by = "sample_name") %>%
    select(species_label, genome_length) %>%
    arrange(species_label)
)