fmu <- '../../data/ciliate_mu_mean_se.csv'
fht <- 'asr_fit_pi4_dist.txt'
fout <- 'species_pop_mu_Ne_maxpi4.csv'

dmu <- read.csv(file = fmu, row.names = 1)
dht <- read.csv(file = fht, sep = ' ', row.names = 1)

## get estimate of min & max mu
pmu <- dmu[grep('Paramecium', rownames(dmu)), ]
rownames(pmu) <- gsub('Paramecium ', '', rownames(pmu))
pmu$maxu <- pmu$mean + 2*pmu$SE
pmu$minu <- pmu$mean - 2*pmu$SE
minu <- min(pmu$minu)
maxu <- max(pmu$maxu)
avgu <- mean(pmu$mean)

dht$sps <- sapply( strsplit(x = rownames(dht), split = '_', fixed = T), function(x)
  paste( head(x,-1), collapse = '_'))
dht$pop <- sapply( strsplit(x = rownames(dht), split = '_', fixed = T), tail, 1)
n <- nrow(dht)

pmu <- pmu[ dht$sps, ]

df <- data.frame( species = dht$sps, population = dht$pop,
                  het = dht$dmax, avgu = pmu$mean, 
                  minu = pmu$minu, maxu = pmu$maxu)

ixs <- is.na(df$avgu)
df[ ixs, 'avgu'] <- avgu
df[ ixs, 'minu'] <- minu
df[ ixs, 'maxu'] <- maxu

## function to get Ne from  mu and H
getNe <- function(u,H,k=4){
  #t <- H/(1-H)
  Ne <- H/(k*u)
  return(Ne)
}

### Ne in aurelia
adf <- df[ grep( "aurelia", df$species), ]
adf$Ne <- apply( adf[ , c('avgu','het')], 1, function(x) getNe(x[1],x[2],k = 2))
adf$minNe <- apply( adf[ , c('maxu','het')], 1, function(x) getNe(x[1],x[2],k = 2))
adf$maxNe <- apply( adf[ , c('minu','het')], 1, function(x) getNe(x[1],x[2],k = 2))

### Ne for others
odf <- df[ grep( "aurelia", df$species, invert = T), ]
odf$Ne <- apply( odf[ , c('avgu','het')], 1, function(x) getNe(x[1],x[2]))
odf$minNe <- apply( odf[ , c('maxu','het')], 1, function(x) getNe(x[1],x[2]))
odf$maxNe <- apply( odf[ , c('minu','het')], 1, function(x) getNe(x[1],x[2]))

df <- as.data.frame( rbind( adf, odf))

df <- df[ order( df$Ne, decreasing = T), ]

## write table
df$avgu <- df$avgu/1e-11
df$minu <- df$minu/1e-11
df$maxu <- df$maxu/1e-11
df$Ne <- df$Ne/1e7
df$minNe <- df$minNe/1e7
df$maxNe <- df$maxNe/1e7

write.csv( x = df, file = fout, quote = F, row.names = F)