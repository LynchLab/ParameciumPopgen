indr <- 'data/pi4_by_dist'
fin <- 'data/pi4_dist_asr_fit.txt'
fout <- 'images/figure_2B.png'

# maximum position used in plotting
maxpos <- 150
# vector of positions used
sitevec <- 1:maxpos

fitfunc <- function(dmin,dmax,c,x=sitevec){
  dmax + (dmin-dmax)*exp(-c*x)
}

# load ASR fit results
dft <- read.csv( file = fin, sep = ' ', row.names = 1)

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
png( fout, width = 6, height = 3, units = 'in', res = 150, pointsize = 8)
par( mar = c(5,5,1,15))
cxl = 1.4
cxx = 1.3
cxp = 0.5

# Set color palette
fcol <- 'species_colors_pts.csv'
dcol <- read.csv(file = fcol,row.names = 1)
cols <- dcol[ gsub('_pop.*', '', sps), 'hex']
cols[is.na(cols)] <- dcol[ 'tetraurelia', 'hex']

# points
ptypes <- dcol[ gsub('_pop.*', '', sps), 'pch']
ptypes[is.na(ptypes)] <- dcol[ 'tetraurelia', 'pch']

# line types
ltypes <- rep( 1, nsp)
ltypes[ grep( 'cladeB', sps)] <- 3

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
mtext(text = 'Distance (bp) from TIS', side = 1, cex = cxl, line = 3)

ytks <- -4:-1
axis( side = 2, at = ytks, labels = as.expression( sapply( ytks, function(x)
  bquote( paste( 10^.(x))))), cex.axis = cxx)
mtext(text = bquote( paste( 'Mean ', pi[4])), side = 2, cex = cxl, line = 3)
box()
## grids
abline( v = xtks, lty = 3, col = 'grey')
abline( h = ytks, lty = 3, col = 'grey')

## add species labels
spnames <- gsub('_pop.+', '', sps)
sshort <- sapply( strsplit( x = sps, split = '', fixed = T), function(x) 
  paste(x[1:4],collapse = ''))
snames <- gsub("_clade", " clade ", spnames)
snames <- gsub("^", "P. ", snames)
snames <- gsub("_.*", "", snames)

legs <- as.data.frame( cbind( snames, cols, ptypes, ltypes))
legs <- legs[ !duplicated(legs), ]

if(F){
labnames <- as.list(legs$snames)
cxs <- grep("clade",legs$snames)
ncxs <- grep("clade",legs$snames,invert = T)
labnames[ncxs] <- sapply( labnames[ncxs], function(x) 
  bquote( paste( italic( .(x)))))
afunc <- function(aname){
  terms <- unlist( strsplit( x = aname, split = ' ', fixed = T))
  bquote( paste( italic( .(terms[1])), ' ', italic( .(terms[2])), ' clade ', .(terms[4])))
}
labnames[cxs] <- sapply( labnames[cxs], afunc)
}

legend( 'topright', 
        legend = legs$snames, col = legs$cols, pch = as.integer(legs$ptypes), 
       lty = as.integer(legs$ltypes), text.font = 3, 
       inset = c(-0.57,-.03), xpd = T, y.intersp = 1.05,
       cex = 1.2, bty = 'n')


dump <- dev.off()