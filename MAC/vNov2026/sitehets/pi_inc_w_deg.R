
fin <- 'species_pop_class_het.csv'
fout <- 'figure_S5.png'

df <- read.csv( file = fin, header = T)

df <- df[ df$population == 'all', ]


## long to wide
df <- df[ , c('species','siteclass','het')]
df <- reshape(data = df, idvar = 'species', timevar = 'siteclass', 
              v.names = 'het', direction = 'wide')

colnames(df) <- gsub( 'het.', '', colnames(df), fixed = T)
rownames(df) <- df$species
########

df <- df[ , c("deg1","deg2","deg3","deg4")]

species <- rownames(df)
nsp <- length(species)

cors <- apply( df, 1, function(x) cor(x = 1:4, y = x, method = 's'))

df <- log10(df)

png( fout, width = 6, height = 4, units = 'in', pointsize = 8, 
     res = 150)
par(mar = c(5,5,2,13))

cxp = 1.4
cxl = 1.3
cxx = 1.2

lytks <- seq(-4,-1,1)

# color
dcol <- read.csv('../../species_colors_pts.csv', row.names = 1)
df$col <- dcol[ rownames(df), 'color']
# points
df$pch <- dcol[ rownames(df), 'pch']

for( i in 1:nsp ){
  if( i == 1)
    plot( 1:4, df[ i, 1:4], ylim = c(-4,-1), col = df$col[i], pch = df$pch[i], cex =cxp, xlab = 'n-fold degeneracy', 
          ylab = bquote(paste(pi[n])), cex.lab = cxl, axes = F, lwd = 1.5)
  else
    points( 1:4, df[ i, 1:4], col = df$col[i], pch = df$pch[i], cex=cxp, lwd = 1.5)
}

# axes
axis( side = 1, at = 1:4, cex.axis = cxx)
axis(side = 2, at = lytks, cex.axis = cxx, labels = as.expression( sapply( lytks, function(x) bquote( paste( 10^.(x))))))
box()

## legend
## add species labels
df$snames <- gsub("_clade", " clade ", rownames(df))
df$snames <- gsub("^", "P. ", df$snames)
df$snames <- gsub("_.*", "", df$snames)

legend( 'topright', 
        legend = df$snames, col = df$col, 
        pch = df$pch, text.font = 3, 
        inset = c(-0.46,-.03), xpd = T, y.intersp = 1.05,
        cex = 1.2, bty = 'n', pt.lwd = 1.5)

dev.off()