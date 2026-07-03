#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

## ----------------------------
## INPUT / OUTPUT
## ----------------------------
FIN <- "~/script/pnps_mean_summary_all_species_no_missing_ratio.csv"


 # if running where file is mounted

OUT_DIR <- "C:/Users/telka/Documents/script"  # change if needed
FOUT    <- file.path(OUT_DIR, "pnps_gene_panel_with_group_boxes.png")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## ----------------------------
## SPECIES COLORS (as provided)
## ----------------------------
col_map <- c(
  burs="#f76666", mult="#9874f4", foki="#ffe599", caud="#e69138",
  sexa="#559738", sonn="#8e7cc3", jenn="#b6d7a8", tetr="#0b5394",
  octa="#a64d79", dode="#b4f4e9", deca="#a2c4c9", tred="#f1c232",
  quad="#ccd9c7", psch="#aaaaaa", quin="#d5ed8a", biau="#ee79ac",
  nov ="#ebc553", tria="#d5a6bd", prim="#6fa8dc", pent="#a95050",
  tbur="#7f5539"
)

to_short <- function(x) {
  x <- tolower(trimws(x))
  if (x == "multi_foki_original") return("mult")
  if (x == "bursaria")            return("burs")
  if (x == "tribursaria")         return("tbur")
  if (x == "sonneborni")          return("sonn")
  if (x == "triaurelia")          return("tria")
  if (x == "primaurelia")         return("prim")
  substr(x, 1, 4)
}

## ----------------------------
## GENE SETS (exact lists; anything else = Unannotated)
## ----------------------------
lgnames <- list(
  ETC = c("nadh1-a","nadh1-b","nadh2","nadh3","nadh4","nadh4L","nadh5",
          "nadh6","nadh7","nadh9","nadh10","cob","cox1","cox2","atp9","yejR"),
  RP  = c("rpL2","rpL6","rpL14","rpL16","rpL32",
          "rpS3","rpS12","rpS13","rpS14","rpS19"),
  YMF = c("ymf56","ymf57","ymf59","ymf61","ymf63","ymf64","ymf65","ymf66",
          "ymf67","ymf68","ymf76","ymf78","ymf79","ymf80","ymf83","ymf84")
)

gene2class <- setNames(
  rep(names(lgnames), lengths(lgnames)),
  unlist(lgnames, use.names = FALSE)
)

## ----------------------------
## LOAD + CLEAN (ignore 0/Inf/NA; force numeric)
## ----------------------------
df <- read_csv(FIN, show_col_types = FALSE) %>%
  mutate(
    Species   = trimws(as.character(Species)),
    Gene      = trimws(as.character(Gene)),
    pN_over_pS = suppressWarnings(as.numeric(pN_over_pS))
  )

df$pN_over_pS[is.infinite(df$pN_over_pS)] <- NA

df <- df %>%
  filter(!is.na(Species), !is.na(Gene), !is.na(pN_over_pS)) %>%
  filter(pN_over_pS > 0)

df <- df %>%
  mutate(
    Class = gene2class[Gene],
    Class = ifelse(is.na(Class), "Unannotated", Class)
  )

## ----------------------------
## SPECIES ORDER + COLORS
## ----------------------------
sps <- sort(unique(df$Species))
nsp <- length(sps)

spshort <- vapply(sps, to_short, character(1))
sp_cols <- unname(col_map[spshort])
sp_cols[is.na(sp_cols)] <- "#555555"

cols <- data.frame(
  sp = sps,
  sp_short = spshort,
  col = sp_cols,
  stringsAsFactors = FALSE
)

## ----------------------------
## GENE ORDER (ETC, RP, YMF, Unannotated; only genes present)
## ----------------------------
genes_present <- sort(unique(df$Gene))

ETC_genes <- intersect(lgnames$ETC, genes_present)
RP_genes  <- intersect(lgnames$RP,  genes_present)
YMF_genes <- intersect(lgnames$YMF, genes_present)

annotated <- unique(c(ETC_genes, RP_genes, YMF_genes))
UN_genes  <- sort(setdiff(genes_present, annotated))

gene_order <- c(ETC_genes, RP_genes, YMF_genes, UN_genes)
df <- df %>% filter(Gene %in% gene_order)

## ----------------------------
## WIDE MATRIX: rows species, cols genes (numeric) -> log10
## ----------------------------
glist <- split(df, df$Gene)
glist <- glist[gene_order[gene_order %in% names(glist)]]
gene_order <- names(glist)

gpnps <- sapply(glist, function(X) {
  v <- X$pN_over_pS[match(sps, X$Species)]
  as.numeric(v)
})
lgpnps <- log10(gpnps)

## ----------------------------
## GROUP BOUNDS (for shading + group boxes)
## ----------------------------
gene_y <- setNames(seq_along(gene_order), gene_order)
ng <- length(gene_order)

grp_genes <- list(
  ETC = intersect(ETC_genes, gene_order),
  RP  = intersect(RP_genes,  gene_order),
  YMF = intersect(YMF_genes, gene_order),
  Unannotated = intersect(UN_genes, gene_order)
)
grp_genes <- grp_genes[sapply(grp_genes, length) > 0]
grp_names <- names(grp_genes)

grp_start <- grp_end <- integer(length(grp_names))
for(i in seq_along(grp_names)){
  ys <- gene_y[grp_genes[[grp_names[i]]]]
  grp_start[i] <- min(ys)
  grp_end[i]   <- max(ys)
}

## ----------------------------
## X LIMITS (log10 scale)
## ----------------------------
xvals <- as.numeric(lgpnps)
xvals <- xvals[is.finite(xvals)]
xmin <- min(-3, floor(min(xvals)))
xmax <- max( 1, ceiling(max(xvals)))

## ----------------------------
## PLOT
## ----------------------------
png(FOUT, width = 14, height = 10, units = "in", res = 150)

# bigger outer box / axes look: thicker lines, larger right margin for labels
op <- par(mar = c(5.5, 8.0, 2.5, 12.5), lwd = 1.6,pin = c(8, 6))

plot(
  NA, NA,
  xlim = c(xmin, xmax + 0.5),
  ylim = c(1, ng),
  xlab = "",
  ylab = "",
  axes = FALSE
)

## Background shading by group block (alternating)
for(i in seq_along(grp_names)){
  rect(
    xleft   = par("usr")[1],
    xright  = par("usr")[2],
    ybottom = grp_start[i] - 0.5,
    ytop    = grp_end[i] + 0.5,
    col     = if (i %% 2 == 1) "grey95" else "white",
    border  = NA
  )
}

# grid
abline(h = 1:ng, lty = 3, col = "grey85")

## ----------------------------
## GROUP BOXPLOTS (bigger + black border + spacing between groups)
## ----------------------------
grp_vecs <- lapply(grp_names, function(g) as.numeric(lgpnps[, grp_genes[[g]], drop = FALSE]))
names(grp_vecs) <- grp_names
bp <- boxplot(x = grp_vecs, plot = FALSE)

gap <- 0.35  # increase for more space between group blocks

for(i in seq_along(grp_names)){
  yb <- grp_start[i] - 0.5 + gap
  yt <- grp_end[i]   + 0.5 - gap
  midy <- (yb + yt)/2
  
  # IQR box (Q1..Q3): larger look + black border
  polygon(
    x = bp$stats[c(2,4,4,2), i],
    y = c(yb, yb, yt, yt),
    col = "grey75",
    border = "black",
    lwd = 1.8
  )
  
  # whiskers (black, thicker)
  segments(bp$stats[1,i], midy, bp$stats[2,i], midy, lwd = 1.8)
  segments(bp$stats[4,i], midy, bp$stats[5,i], midy, lwd = 1.8)
  
  # caps
  segments(bp$stats[1,i], midy - 0.25, bp$stats[1,i], midy + 0.25, lwd = 1.8)
  segments(bp$stats[5,i], midy - 0.25, bp$stats[5,i], midy + 0.25, lwd = 1.8)
  
  # median (bold vertical line spanning the group band)
  segments(bp$stats[3,i], yb, bp$stats[3,i], yt, lwd = 2.6)
}

## ----------------------------
## POINTS PER GENE ROW (species-colored)
## ----------------------------
for(j in seq_len(ng)){
  gene <- gene_order[j]
  pts <- lgpnps[, gene]
  ok  <- is.finite(pts)
  if(any(ok)){
    points(
      x = pts[ok],
      y = rep(j, sum(ok)),
      pch = 21,
      bg  = cols$col[ok],
      col = "grey10",
      cex = 1.45
    )
  }
}

## Gene labels (left)
text(
  x = par("usr")[1] + 0.09*(xmax-xmin),
  y = 1:ng,
  labels = gene_order,
  adj = 1,
  xpd = TRUE,
  cex = 0.82,
  font = 3
)

## X-axis ticks as 10^x
xtks <- seq(xmin, xmax, by = 1)
axis(
  side = 1,
  at = xtks,
  labels = as.expression(sapply(xtks, function(x) bquote(10^.(x)))),
  lwd = 1.6
)
mtext(text = expression(pi[N] / pi[S]), side = 1, line = 3.2, cex = 1.3)

## Group labels on the right
pretty_grp <- c(
  ETC = "Electron transport\nchain",
  RP  = "Ribosomal\nproteins",
  YMF = "Ciliate-specific\ngenes",
  Unannotated = "Unannotated\ngenes"
)

for(i in seq_along(grp_names)){
  ymid <- (grp_start[i] + grp_end[i]) / 2
  text(
    x = par("usr")[2] + 0.12,
    y = ymid,
    labels = pretty_grp[grp_names[i]],
    adj = 0,
    xpd = TRUE,
    cex = 1.45
  )
}

## Legend
legend(
  "topleft",
  inset = c(-0.27, 0),
  legend = cols$sp_short,
  fill = cols$col,
  border = "grey20",
  cex = 1.05,
  xpd = TRUE,
  bty = "n"
)

## Outer box thicker / bigger look
box(lwd = 6)

par(op)
dev.off()

cat("✅ Wrote: ", FOUT, "\n", sep = "")
