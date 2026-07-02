ffn <- 'data/gene_cat_class.tsv'
fni <- 'data/gene_sp_ni_pairwise.csv'
fout <- 'images/figure_4D.png'

# minimum no. of entries in a category
mindata <- 15

dfn <- read.csv(file = ffn, sep = '\t')
dni <- read.csv(file = fni)

# remove category - others
dfn <- dfn[ !dfn$category %in% c("others","development"), ]

# take log
dni$NI <- log10(dni$NI)

# list of NI for each gene
ls_gene_ni <- split( dni$NI, dni$gene)

# limit data to genes common to both datasets
genes <- intersect( dfn$gene, names(ls_gene_ni))
ls_gene_ni <- ls_gene_ni[ genes]
dfn <- dfn[ dfn$gene %in% genes, ]

# group genes by category
cat_genes <- split( dfn$gene, dfn$category)


cat_ni <- lapply( cat_genes, function(x) unlist(ls_gene_ni[x]))
cat_ni <- cat_ni[ sapply( cat_ni, length) >= mindata]
ncat <- length(cat_ni)

med_ni <- sort( sapply( cat_ni, median))
cat_ni <- cat_ni[ names(med_ni)]

cat_ngenes <- sapply( cat_ni, function(x) length( unique( sapply( strsplit( x = names(x), split = '', fixed = T), function(y) paste(y[1:8],collapse = '')))))

cat_ng <- sapply( cat_ni, length)

wtps <- rep( NA, ncat)
names(wtps) <- names(cat_ni)

## tests
for( i in 1:ncat){
  wtps[i] <- wilcox.test( x = cat_ni[[i]], mu = 0, alternative = 'g')$p.value
}

wtps.adj <- p.adjust(p = wtps, method = 'holm')

nsig <- sum( wtps.adj < 0.05)
cat(sprintf("# Categories with NI significantly greater than 1 = %i / %i\n",nsig,ncat))

## plotting
png( fout, width = 4.3, height = 3.4, units = 'in', res = 150, pointsize = 8)
par( mar = c(7.5,4.3,1.5,1))

cxl = 1.4
cxx = 1.3
cxp = 0.8
ptype = 1

ylims <- c(-1.5,1.5)
lytks <- -2:2

bp <- boxplot( cat_ni, axes = F, ylim = range(ylims),
               ylab = "Neutrality Index", boxwex = 0.7, 
               cex.lab = cxl)

abline( h = 0, lty = 2)

tick <- seq_along(bp$names)
op <- par(family = "Helvetica Narrow")
axis(1, at = tick, labels = F)
text( tick, par("usr")[3] - 0.2, gsub("_", " ", bp$names), srt = 45, xpd = T, 
      cex = 1.1, adj = 1)

axis( side = 2, at = lytks, cex.axis = cxx,
      labels = as.expression( sapply( lytks, function(x) 
        bquote( paste( 10^.(x))))))

#axis(side = 3, at = tick, labels = F)
#text( tick, par("usr")[4] + 0.2, cat_ng, xpd = T, cex = 0.7)
box()

dump <- dev.off()