fin <- 'species_pop_class_het.csv'
fout <- 'figure_3A.svg'

df <- read.csv( file = fin)

## long to wide
df <- df[ , c('species','population','siteclass','het')]
df <- reshape(data = df, idvar = c('species','population'), timevar = 'siteclass', 
              v.names = 'het', direction = 'wide')

colnames(df) <- gsub( 'het.', '', colnames(df), fixed = T)
########

df$p1p4 <- df$deg1/df$deg4
df$lp1p4 <- log10(df$p1p4)

svg( fout, width = 2.5, height = 2.5)
par( mar = c(3.3,3.3,1.5,0.8))

cxl = 1.0
cxx = 1.0
cxp = 1.2
cxt = 0.9

ptype = 16

tks <- c( 0.05, 0.1, 0.25)
ltks <- log10(tks)

cols <- read.csv('../../species_colors_pts.csv', row.names = 1)
df$cols <- cols[ df$species, 'color']
# points
df$pch <- cols[ df$species, 'pch']

# species names
df$sp <- gsub("_.*", "", df$species)
sclade <- sapply( strsplit( x = df$species, split = '_clade', fixed = T), '[', 2)
snames <- apply( cbind( df$sp, sclade), 1, paste, collapse='-')
df$snames <- gsub(pattern = '-NA$', replacement = '', x = snames)

# separate data for the sets of all isolates pooled
alldf <- df[ df$population == 'all', c("species", "lp1p4","snames","cols","pch")]
df <- df[ df$population != 'all', ]
df$all <- alldf$lp1p4[ match( df$species, alldf$species)]

df <- df[ !is.infinite(df$lp1p4), ]

plot( lp1p4 ~ all, data = df, 
      xlab = "", 
      ylab = "", 
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
#mtext( text = bquote( paste( italic(r), ' = ', .(r))), 
#       side = 3, adj = 0, line = -1, cex = cxx)

## axes
axis( side = 1, at = ltks, labels = tks, cex.axis = cxx)
axis( side = 2, at = ltks, labels = tks, cex.axis = cxx)
box()
mtext( text = bquote( paste( 'Pooled isolates ', pi[0]/pi[4])),
       side = 1, line = 2.4, cex = cxl)
mtext( text = bquote( paste( 'Population ', pi[0]/pi[4])), side = 2, 
       line = 2, cex = cxl)


legs <- alldf[ , c('snames', 'cols','pch')]
legs <- legs[ order(legs$snames), ]

dump <- dev.off()