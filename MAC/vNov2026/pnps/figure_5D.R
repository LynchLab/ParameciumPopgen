ffn <- 'gene_cat_class.tsv'
fni <- 'gene_sp_ni_pairwise.csv'
fout <- 'figure_5D.svg'

# minimum no. of entries in a category
mindata <- 15

# maximum dS for a comparison
maxds <- 0.8

dfn <- read.csv(file = ffn, sep = '\t')
dni <- read.csv(file = fni)

# apply max dS threshold
dni <- dni[ dni$dS <= maxds, ]

# species with valid NI 
sps <- unique( c(dni$focal_species,dni$closest_species))
cat( sprintf("Species with valid NI = %s\n", paste(sps,collapse = ', ')))

# take log
dni$NI <- log10(dni$NI)

# list of NI for each gene
ls_gene_ni <- split( dni, dni$gene)

posg <- names(ls_gene_ni)[sapply( ls_gene_ni, function(X) all(X$NI<0))]
writeLines(text = posg, con = "genes_all_ni_lt1.txt")



# remove category - others
dfn <- dfn[ !dfn$category %in% c("uncategorized"), ]





# limit data to genes common to both datasets
genes <- intersect( dfn$gene, names(ls_gene_ni))
cat( sprintf("# genes with valid NI = %s\n", length(genes)))

ls_gene_ni <- ls_gene_ni[ genes]
dfn <- dfn[ dfn$gene %in% genes, ]

# group genes by category
cat_genes <- split( dfn$gene, dfn$category)
cat(sprintf("# Gene categories in the data = %i\n",length(cat_genes)))

cat_ni <- lapply( cat_genes, function(x) 
  do.call( rbind, ls_gene_ni[x]))
cat_ni <- cat_ni[ sapply( cat_ni, nrow) >= mindata]
ncat <- length(cat_ni)
cat(sprintf("# Gene categories with at least %i comparisons = %i\n", mindata, ncat))

# genes within these categories
genes <- unlist( sapply( cat_ni, function(x) unique(x$gene)))
cat( sprintf("# genes within these categories = %s\n", length(genes)))

# sort categories by median
med_ni <- sort( sapply( cat_ni, function(X) median(X$NI)))
cat_ni <- cat_ni[ names(med_ni)]

cat_ngenes <- sapply( cat_ni, function(x) length( unique(x$gene)))

cat_ng <- sapply( cat_ni, length)

## Wilcoxon tests
wtps <- rep( NA, ncat)
names(wtps) <- names(cat_ni)
for( i in 1:ncat){
  wtps[i] <- wilcox.test( x = cat_ni[[i]][,'NI'], mu = 0, alternative = 'g')$p.value
}

wtps.adj <- p.adjust(p = wtps, method = 'holm')

nsig <- sum( wtps.adj < 0.05)
cat(sprintf("# Categories with NI significantly greater than 1 = %i\n",nsig))

## plotting
cat_ni_only <- sapply( cat_ni, function(X) X[,'NI'])

svg( fout, width = 4, height = 3)
par( mar = c(6,3.6,1.5,1))

cxl = 1.0
cxx = 1.0
cxt = 0.75
cxp = 0.8

ptype = 1

ylims <- c(-1.5,1.5)
lytks <- -2:2

acol <- 'darkslategrey'
tcol <- apply( col2rgb(col = acol), 2, function(x) 
  rgb( x[1], x[2], x[3], 100, maxColorValue = 255))


bp <- boxplot( cat_ni_only, axes = F, ylim = range(ylims),
               ylab = "", boxwex = 0.6, cex = cxp, col = tcol,
               cex.lab = cxl, border = acol, lwd = 0.8, outline = T)

abline( h = 0, lty = 2)

mtext( text = "Neutrality Index", 
       side = 2, line = 2.5, cex = cxl)


## labels
tick <- seq_along(bp$names)
#op <- par(family = "Helvetica Narrow")
axis(1, at = tick, labels = F)
text( tick, par("usr")[3] - 0.3, gsub("_", " ", bp$names), srt = 45, xpd = T, 
      cex = cxt, adj = 1)

axis( side = 2, at = lytks, cex.axis = cxx, las = 2,
      labels = as.expression( sapply( lytks, function(x) 
        bquote( paste( 10^.(x))))))


#axis(side = 3, at = tick, labels = F)
#text( tick, par("usr")[4] + 0.2, cat_ng, xpd = T, cex = 0.7)

box()

dump <- dev.off()