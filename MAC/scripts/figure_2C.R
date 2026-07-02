  fne <- 'data/species_pop_mu_Ne_maxpi4_asr.csv'
fpi <- 'data/species_pop_class_het_gfe.csv'
fout <- 'images/figure_2C.png'

dne <- read.csv(file = fne)
dne$Ne <- dne$Ne*1e07
dne$lNe <- log10(dne$Ne)

df <- read.csv( file = fpi)
df <- df[ df$population != 'all', ]
df$lp1p4 <- log10(df$p1p4)

rownames(df) <- apply( df[,c('species','population')], 1, paste, collapse='_')
rownames(dne) <- apply( dne[,c('species','population')], 1, paste, collapse='_')

df$lNe <- dne[ match( rownames(df), rownames(dne)), 'lNe']

df <- df[ order(df$lNe), ]

## check correlation between Ne ( from asr pi4 max) and theta 4 (angsd)
cor.test( df$lNe, df$deg4, method = 's')

png( fout, width = 3, height = 3, units = 'in', res = 150, pointsize = 8)
par( mar = c(4.5,4.5,2,2))

cxl = 1.4
cxx = 1.3
cxp = 1.5
ptype = 16

# axis ticks
xtks <- 6:9
ytks <- c(0.02, 0.05, 0.2, 0.5)
lytks <- log10(ytks)

cols <- read.csv('species_colors_pts.csv', row.names = 1)

# color
df$cols <- cols[ df$species, 'color']
# points
df$pch <- cols[ df$species, 'pch']

# scatterplot
plot( lp1p4 ~ lNe, data = df, xlab = expression(italic(N)[e]),
      ylab = expression(pi[0]/pi[4]), col = df$cols, pch = df$pch, cex = cxp,
      cex.lab = cxl, axes = F, xlim = range(xtks), ylim = range(lytks))

# linear regression
lob <- lm( lp1p4 ~ lNe, data = df)
abline( lob, lty = 2)
lob.sm <- summary(lob)
b <- round( lob.sm$coefficients[2,1], 2)
d <- round( lob.sm$coefficients[2,2], 2)
mtext( text = bquote( paste( beta, ' = ', .(b), ' \u00B1 ', .(d))), side = 3, adj = 1, cex = cxx)

# spearman correlation
cob <- cor.test( df$lNe, df$lp1p4, method = 's')
r <- round( cob$estimate, 2)
mtext( text = bquote( paste( italic(r), ' = ', .(r))), side = 3, adj = 0, cex = cxx)

# axes
axis(side = 1, at = xtks, labels = as.expression( 
  sapply( xtks, function(x) bquote( paste( 10^.(x))))), cex.axis = cxx)
axis(side = 2, at = lytks, labels = ytks, cex.axis = cxx)
box()

dump <- dev.off()