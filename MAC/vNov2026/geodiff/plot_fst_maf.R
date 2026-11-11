
indr <- "PSA"
fout <- "fst_stats.csv"
fplot <- "figure_S4.png"

fls <- list.files(indr)

nf <- length(fls)
## extract pair names
pnames <- gsub( ".txt", "", fls)
pnames <- gsub( "_", " ", pnames)
pnames <- gsub(" clade", "-", pnames)
pnames <- gsub("pop", "/ ", pnames)

wsize <- 0.025
fstart <- seq( 0, 0.5-wsize, wsize)
fbins <- cbind( fstart, fstart + wsize)
nb <- nrow(fbins)

mfstmat <- matrix( NA, nb, nf+2)
mfstmat[,1:2] <- fbins

colnames(mfstmat) <- c("min","max", pnames)

qts <- matrix( NA, nf, 5, 
               dimnames = list( pnames, c('25%','50%','75%','mean','se')))

for( i in 1:nf){
  df <- read.csv( file = file.path( indr, fls[i]), sep = '\t')
  df <- df[ !is.na(df$theta), ]
  df <- df[ df$theta >= 0, ]
  qts[i,1:3] <- quantile(x = df$theta, probs = c(0.25,0.5,0.75))
  qts[i,4] <- mean(df$theta)
  qts[i,5] <- sd(df$theta)/sqrt(nrow(df)-1)
  mfstmat[,i+2] <- apply( fbins, 1, function(x) mean( df[ df$MAF > x[1] & df$MAF <= x[2], 'theta']))
}

write.csv( x = qts, file = fout)

## mean over last windows (plateau)
platfst <- colMeans(mfstmat[ mfstmat[,1]>=0.4,-(1:2)], na.rm = T)
platfst <- platfst[ order(platfst, decreasing = T)]


## plotting
png( fplot, width = 6, height = 6, units = 'in', res = 150, pointsize = 10)
par( mfrow = c(3,3), mar = c(0,0,0,0), oma = c(5,5.5,2,2))

cxx = 1.2
cxl = 1.2

for( i in 1:nf){
  # mean Fst over frequency bins
  plot( mfstmat[,2], mfstmat[,2+i], 
        xlab = '', ylab = '',
        axes = F, xlim = c(0,0.5), ylim = c(0,1))
  box()
  # quantiles
  abline( h = qts[i,2], lty = 2)
  abline( h = qts[i,c(1,3)], lty = 3)
  
  # title
  mtext(text = pnames[i], side = 3, line = -2.5, adj = 0.9)
  
  # x-axis
  if(i %in% 7:9){
    axis(side = 1, at = seq(0.1,0.5,0.1), cex.axis = cxx)
  }
  # y-axis
  if(i %in% c(1,4,7)){
    axis(side = 2, at = seq(0.2,1,0.2), cex.axis = cxx)
  }
}

# axis labels
mtext(text = "Minor Allele Frequency", side = 1, line = 3.5, outer = T, cex = cxl)
mtext(text = expression(italic(F)[ST]), side = 2, line = 3, outer = T, cex = cxl)

dev.off()