## Plot p0/p4 against Ne (estimated from p4 max) using different pops for each
fne <- 'species_pop_mu_Ne_maxpi4.csv'
fpi <- 'species_pop_class_het.csv'
fout <- 'figure_S6.png'

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

spdfs <- split( df, df$species)
spdfs <- spdfs[ sapply( spdfs, nrow) >= 2]

nsp <- length(spdfs)

out <- NULL

for( i in 1:nsp){
  X <- spdfs[[i]]
  X <- X[ order( X$lNe, decreasing = T), ]
  sdf <- data.frame( species = X$species[-1], Ne = X$lNe[1], pnps = X$lp1p4[-1])
  out <- rbind( out, sdf)
}

df <- out
df <- df[ !is.na(df$Ne), ]
df <- df[ order(df$Ne), ]

## plotting
png( fout, width = 4.5, height = 3, units = 'in', pointsize = 8, 
     res = 150)
par(mar = c(4.2,4.5,2,10.5))

cxl = 1.3
cxx = 1.2
cxp = 1.5

# axis ticks
xtks <- 6:9
ytks <- c(0.02, 0.05, 0.2, 0.5)
lytks <- log10(ytks)

# Set color palette
fcol <- '../../species_colors_pts.csv'
cols <- read.csv( file = fcol, row.names = 1)
df$cols <- cols[ df$species, 'color']
# points
df$pch <- cols[ df$species, 'pch']


# scatterplot
plot( pnps ~ Ne, data = df, xlab = expression(italic(N)[e]),
      ylab = expression(pi[0]/pi[4]), col = df$cols, pch = df$pch, 
      cex = cxp, lwd = 1.5, cex.lab = cxl, 
      axes = F, xlim = range(xtks), ylim = range(lytks))

# regression
lob <- lm( pnps ~ Ne, data = df)
abline( lob, lty = 2)
lob.sm <- summary(lob)
b <- round( lob.sm$coefficients[2,1], 2)
d <- round( lob.sm$coefficients[2,2], 2)
mtext( text = bquote( paste( beta, ' = ', .(b), ' \u00B1 ', .(d))), side = 3, adj = 1, cex = cxx)

# correlation
cob <- cor.test(df$Ne,df$pnps,method='s')
r <- round( cob$estimate, 2)
p <- round( cob$p.value, 4)
mtext( text = bquote( paste( italic(r), ' = ', .(r), ', ', italic(P), ' = ', .(p))), side = 3, adj = 0, cex = cxx)


# axes
axis(side = 1, at = xtks, labels = as.expression( 
  sapply( xtks, function(x) bquote( paste( 10^.(x))))), cex.axis = cxx)
axis(side = 2, at = lytks, labels = ytks, cex.axis = cxx)
box()

## legend
## add species labels
df$snames <- gsub("_clade", " clade ", df$species)
#df$snames <- gsub("^", "P. ", df$snames)
df$snames <- gsub("_.*", "", df$snames)

legs <- df[ , c( "snames", "cols", "pch")]
legs <- unique(legs)

legend( 'topright', 
        legend = legs$snames, col = legs$col, 
        pch = legs$pch, text.font = 3, 
        inset = c(-0.56,-.03), xpd = T, y.intersp = 1.05,
        cex = 1.2, bty = 'n', pt.lwd = 1.5)

dump <- dev.off()