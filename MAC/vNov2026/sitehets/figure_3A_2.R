fin <- 'species_pop_class_het.csv'
fiso <- '../geodiff/isolate_location_coords.csv'
fout <- 'figure_3A.svg'

df <- read.csv( file = fin)

## long to wide
df <- df[ , c('species','population','siteclass','het')]
df <- reshape(data = df, idvar = c('species','population'), timevar = 'siteclass', 
              v.names = 'het', direction = 'wide')

colnames(df) <- gsub( 'het.', '', colnames(df), fixed = T)
rownames(df) <- apply( df[,c("species","population")], 
                                    1, paste, collapse="_")
########

## isolate list ##############
diso <- read.csv( file = fiso)
# no. of isolates by species
sp_niso <- table(diso$species)
# no. of isolates by population in each species
diso <- diso[ diso$population != "", ]
diso$spop <- apply( diso[,c("species","population")], 1, paste, collapse="_pop")
pop_niso <- table(diso$spop)
# keep populations with at least 4 isolates
pop_niso <- pop_niso[ pop_niso >= 4]
# attach number of isolates to main data
df$n <- pop_niso[ match( rownames(df), names(pop_niso))]
df[ grep("all",rownames(df)), 'n'] <-  sp_niso[df$species[grep("all",rownames(df))]]
# no. of pairs
df$np <- choose( df$n, 2)

## attach species data to each population
df <- df[ c("species","population","deg1","deg4","np")]

# separate data for the sets of all isolates pooled
alldf <- df[ df$population == 'all', c("species", "deg1","deg4", "np")]
df <- df[ df$population != 'all', ]

ixs <- match( df$species, alldf$species)
df <- as.data.frame( cbind( df, alldf[ ixs, c("deg1","deg4", "np")]))
colnames(df) <- c("species","population","popdeg1","popdeg4","popnp",
                  "alldeg1","alldeg4","allnp")

# pis excluding pop
df$pix1 <- (df$alldeg1*df$allnp - df$popdeg1*df$popnp)/(df$allnp-df$popnp)
df$pix4 <- (df$alldeg4*df$allnp - df$popdeg4*df$popnp)/(df$allnp-df$popnp)

## piN/piS #############
df$pop_p1p4 <- df$popdeg1/df$popdeg4
df$exl_p1p4 <- df$pix1/df$pix4

## logs
df$lpopp1p4 <- log10(df$pop_p1p4)
df$lexlp1p4 <- log10(df$exl_p1p4)

df <- df[ !is.infinite(df$lpopp1p4), ]

### plotting
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

x <- as.numeric( df$lexlp1p4)
y <- as.numeric( df$lpopp1p4)

plot( x, y,
      xlab = "",ylab = "", 
      col = df$cols, pch = df$pch, cex = cxp,
      cex.lab = cxl, axes = F, xlim = range(ltks), ylim = range(ltks))

lob <- lm( lpopp1p4 ~ lexlp1p4, data = df)
abline( lob, lty = 2)
abline( 0, 1, lty = 3)

lob.sm <- summary(lob)
b <- round( lob.sm$coefficients[2,1], 2)
d <- round( lob.sm$coefficients[2,2], 2)
mtext( text = bquote( paste( beta, ' = ', .(b), ' \u00B1 ', .(d))), side = 3, adj = 1, cex = cxx)

cob <- suppressWarnings( cor.test( df$pop_p1p4, df$exl_p1p4, method = 's'))
r <- round( cob$estimate, 2)

## axes
axis( side = 1, at = ltks, labels = tks, cex.axis = cxx)
axis( side = 2, at = ltks, labels = tks, cex.axis = cxx)
box()
mtext( text = bquote( paste( 'Pooled isolates ', pi[0]/pi[4])),
       side = 1, line = 2.4, cex = cxl)
mtext( text = bquote( paste( 'Population ', pi[0]/pi[4])), side = 2, 
       line = 2, cex = cxl)

dump <- dev.off()