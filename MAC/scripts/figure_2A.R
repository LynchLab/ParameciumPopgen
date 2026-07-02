fin <- 'data/species_pop_class_het_gfe.csv'
fout <- 'images/figure_2A.png'

df <- read.csv( file = fin)
df$lp1p4 <- log10(df$p1p4)

png( fout, width = 3, height = 3, units = 'in', res = 150, pointsize = 8)
par( mar = c(4.5,4.5,2,2))

cxl = 1.4
cxx = 1.3
cxp = 1.5
ptype = 16

tks <- c( 0.05, 0.1, 0.25)
ltks <- log10(tks)

cols <- read.csv('species_colors_pts.csv', row.names = 1)
df$cols <- cols[ df$species, 'color']
# points
df$pch <- cols[ df$species, 'pch']

# species names
sclade <- sapply( strsplit( x = df$species, split = '_clade', fixed = T), '[', 2)
snames <- apply( cbind( df$sp, sclade), 1, paste, collapse='-')
df$snames <- gsub(pattern = '-NA$', replacement = '', x = snames)

# separate data for the sets of all isolates pooled
alldf <- df[ df$population == 'all', c("species", "lp1p4","snames","cols","pch")]
df <- df[ df$population != 'all', ]
df$all <- alldf$lp1p4[ match( df$species, alldf$species)]

plot( lp1p4 ~ all, data = df, xlab = bquote( paste( 'Pooled isolates ', pi[0]/pi[4])), 
      ylab = bquote( paste( 'Population ', pi[0]/pi[4])), 
      col = df$cols, pch = df$pch, cex = cxp,
      cex.lab = cxl, axes = F, xlim = range(ltks), ylim = range(ltks))

lob <- lm( lp1p4 ~ all, data = df)
abline( lob, lty = 2)
abline( 0, 1, lty = 3)

lob.sm <- summary(lob)
b <- round( lob.sm$coefficients[2,1], 2)
d <- round( lob.sm$coefficients[2,2], 2)
mtext( text = bquote( paste( beta, ' = ', .(b), ' \u00B1 ', .(d))), side = 3, adj = 1, cex = cxx)

cob <- suppressWarnings( cor.test( df$all, df$lp1p4, method = 's'))
r <- round( cob$estimate, 2)
mtext( text = bquote( paste( italic(r), ' = ', .(r))), side = 3, adj = 0, cex = cxx)

axis( side = 1, at = ltks, labels = tks, cex.axis = cxx)
axis( side = 2, at = ltks, labels = tks, cex.axis = cxx)
box()

legs <- alldf[ , c('snames', 'cols','pch')]
legs <- legs[ order(legs$snames), ]

dump <- dev.off()