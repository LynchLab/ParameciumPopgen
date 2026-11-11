indr <- 'pi4_by_dist'
fin <- 'asr_fit_pi4_dist.txt'
fout <- 'figure_3B.svg'

# maximum position used in plotting
maxpos <- 150
# vector of positions used
sitevec <- 1:maxpos

fitfunc <- function(dmin,dmax,c,x=sitevec){
  dmax + (dmin-dmax)*exp(-c*x)
}

# load ASR fit results
dft <- read.csv( file = fin, sep = ' ', row.names = 1)
dft <- dft[ !grepl("all", rownames(dft)), ]

# order by decreasing value of dmax
dft <- dft[ order(dft$dmax, decreasing = T), ]
sps <- rownames(dft)
nsp <- length(sps)

## read data, get log mean div of positions up to maxpos
dfs <- sapply( sps, function(X) NULL)
for(i in 1:nsp){
    fpath <- file.path( indr, paste0( sps[i], '.csv'))
    df <- read.csv( fpath, header = T)
    df$ldiv <- log10(df$mdiv)
    df <- df[ df$site %in% sitevec, ]
    dfs[[i]] <- df
}


## plotting
svg( fout, width = 5, height = 2.7)
par( mar = c(3.3,3.5,1,8.5))

cxl = 1.0
cxx = 1.0
cxp = 0.3
cxt = 0.9

# Set color palette
fcol <- '../../species_colors_pts.csv'
dcol <- read.csv(file = fcol,row.names = 1)
cols <- dcol[ gsub('_pop.*', '', sps), 'hex']

# points
ptypes <- dcol[ gsub('_pop.*', '', sps), 'pch']

# line types
ltypes <- rep( 1, nsp)
#ltypes[ grep( 'cladeB', sps)] <- 3

# for each speciees
for(i in 1:nsp){
  # plot log mean div against sites
  if(i == 1){
    plot( ldiv ~ site, data = dfs[[i]], 
          col = cols[i], pch = ptypes[i], cex = cxp, lwd = 0.5, 
          xlim = c(0,maxpos), ylim = c(-3.5,-1), axes = F, xlab = '', ylab = '')
  }else{
    points( ldiv ~ site, data = dfs[[i]], col = cols[i], pch = ptypes[i], cex = cxp, lwd = 0.5)
  }
  # show ASR fit
  predy <- log10(fitfunc(dmin = dft[sps[i],'dmin'], dmax = dft[sps[i],'dmax'],
                   c = dft[ sps[i], 'c']))
  #abline( h = median(dfs[[i]]$ldiv), lty = 1, col = cols[i])
  lines( sitevec, predy, col = cols[i], lwd = 1, lty = ltypes[i])
}
## axes
xtks <- seq(0,maxpos,by=50)
axis( side = 1, at = xtks, cex.axis = cxx)
mtext(text = 'Distance (bp) from TIS', side = 1, cex = cxl, line = 2)

ytks <- -4:-1
axis( side = 2, at = ytks, labels = as.expression( sapply( ytks, function(x)
  bquote( paste( 10^.(x))))), cex.axis = cxx, las = 2)
mtext(text = bquote( paste( 'Mean ', pi[4])), side = 2, cex = cxl, line = 2.4)
box()

## add species labels
spnames <- gsub('_pop.+', '', sps)
sshort <- sapply( strsplit( x = sps, split = '', fixed = T), function(x) 
  paste(x[1:4],collapse = ''))
snames <- gsub("_clade", " clade ", spnames)
snames <- gsub("_.*", "", snames)

legs <- as.data.frame( cbind( snames, cols, ptypes, ltypes))
legs <- legs[ !duplicated(legs), ]

legend( 'topright', 
        legend = legs$snames, col = legs$cols, pch = as.integer(legs$ptypes), 
       lty = as.integer(legs$ltypes), text.font = 3, 
       inset = c(-0.67,-0.04), xpd = T, 
       x.intersp = 0.8, y.intersp = 1,
       cex = 0.85, bty = 'n')


dump <- dev.off()