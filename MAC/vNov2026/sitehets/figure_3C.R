fne <- 'species_pop_mu_Ne_maxpi4.csv'
fpi <- 'species_pop_class_het.csv'
fout <- 'figure_3C.svg'

dne <- read.csv(file = fne)
dne$Ne <- dne$Ne*1e07
dne$lNe <- log10(dne$Ne)

df <- read.csv( file = fpi)
df <- df[ df$population != 'all', ]

## long to wide
df <- df[ , c('species','population','siteclass','het')]
df <- reshape(data = df, idvar = c('species','population'), timevar = 'siteclass', 
              v.names = 'het', direction = 'wide')

colnames(df) <- gsub( 'het.', '', colnames(df), fixed = T)
########

df$lp1p4 <- log10(df$deg1/df$deg4)

rownames(df) <- apply( df[,c('species','population')], 1, paste, collapse='_')
rownames(dne) <- apply( dne[,c('species','population')], 1, paste, collapse='_')

df$lNe <- dne[ match( rownames(df), rownames(dne)), 'lNe']

df <- df[ order(df$lNe), ]

## check correlation between Ne ( from asr pi4 max) and theta 4 (angsd)
cor.test( df$lNe, df$deg4, method = 's')

svg( fout, width = 2.5, height = 2.5)
par( mar = c(3.3,3.3,1.5,0.8))

cxl = 1.0
cxx = 1.0
cxp = 1.2
cxt = 0.9

ptype = 16

# axis ticks
xtks <- 6:9
ytks <- c(0.05, 0.1, 0.25)
lytks <- log10(ytks)

cols <- read.csv('../../species_colors_pts.csv', row.names = 1)

# color
df$cols <- cols[ df$species, 'color']
# points
df$pch <- cols[ df$species, 'pch']

# scatterplot
plot( lp1p4 ~ lNe, data = df, xlab = "",
      ylab = "", col = df$cols, pch = df$pch, cex = cxp,
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
#mtext( text = bquote( paste( italic(r), ' = ', .(r))), side = 3, adj = 0, cex = cxx)

# axes
axis(side = 1, at = xtks, labels = as.expression( 
  sapply( xtks, function(x) bquote( paste( 10^.(x))))), cex.axis = cxx)
axis(side = 2, at = lytks, labels = ytks, cex.axis = cxx)

mtext( text = expression( italic(N)[e]),
       side = 1, line = 2.4, cex = cxl)
mtext( text = expression(pi[0]/pi[4]), side = 2, 
       line = 2, cex = cxl)




box()

dump <- dev.off()