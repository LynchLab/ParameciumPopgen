### Test population differentiation by permuting location identities
dgen <- 'gendist_for_samples_w_coordinates'
floc <- 'isolate_location_coords.csv'

fout <- 'popdiff_perms.csv'

maxnperm <- 10000

tstart <- Sys.time()

sps <- gsub( ".csv", "", list.files(dgen))
nsp <- length(sps)

## extract location ids for each isolate #############
dloc <- read.csv(file = floc)

dloc <- dloc[ !is.na(dloc$latitude), ]

isolocs <- gsub(pattern = " ", replacement = "", x = dloc$source)
names(isolocs) <- dloc$isolate

### Function to get mean distance among isolates at a location
mdist_in_pop <- function(isolist,distmat){
  # number of locations
  nloc <- length(isolist)
  # initialize vector of mean distances
  mds <- rep( NA, nloc)
  # for each location
  for(j in 1:nloc){
    # vector of isolates
    iso <- isolist[[j]]
    # their pairwise distance matrix
    lgen <- distmat[ iso, iso]
    # its mean
    mds[j] <- mean( lgen[ lower.tri(lgen)])
  }
  # mean of means
  mean(mds)
}

### Function to permutate location labels
## while keeping sample size constant for each population
permute_loc <- function(isolist){
  # number of isolates at each location
  nisos <- sapply( isolist, length)
  # location names
  ulocs <- names(isolist)
  # vector of locations with each repeated up to its sample size
  rlocs <- rep( ulocs, nisos)
  # shuffle isolate names and attach to location names
  risos <- sample( unlist(isolist), replace = F)
  split( risos, rlocs)
}

### Function to get mean distances over multiple permutations
get_rdist <- function(isolist,distmat,nperm=1000){
  replicate( n = nperm, expr = mdist_in_pop(isolist = permute_loc(isolist), 
                                            distmat = distmat))
}

### get no. of possible permutations
get_nperm <- function(isolist){
  loc_niso <- sapply( isolist, length)
  niso <- sum(loc_niso)
  exp( lfactorial(niso) - sum(lfactorial(loc_niso)))
}


out <- matrix( NA, nsp, 6, dimnames = list( sps, c( "nloc", "obs_mean", "nperm", "cileft","ciright",
                                                   "sig")))

for(i in 1:nsp){
  cat( sprintf( "%s: species = %s\n",Sys.time(),sps[i]))
  # genetic distance matrix
  fgen <- file.path( dgen, paste0(sps[i],".csv"))
  genmat <- read.csv( fgen, row.names = 1, check.names = F)
  
  # isolates
  isolates <- rownames(genmat)
  locids <- isolocs[ isolates]
  # clear all spaces from location ids
  locids <- gsub(" +", "", locids)
  loc_iso <- split( isolates, locids)
  # remove singletons
  loc_iso <- loc_iso[ sapply( loc_iso, length)>1]
  nloc <- length(loc_iso)
  cat( sprintf("\t # locations w/ multiple isolates = %i\n",nloc))
  
  if( nloc < 2){
    out[ i, 1] <- nloc
    next
  }
  
  obsm <- mdist_in_pop(isolist = loc_iso, distmat = genmat)
  
  ## permutations
  totnp <- get_nperm(isolist = loc_iso)
  cat( sprintf("\t # possible permutations = %.0f\n",totnp))
  nperms <- ifelse( totnp < maxnperm, totnp, maxnperm)
  
  rdrb <- get_rdist(isolist = loc_iso, distmat = genmat, 
                    nperm = nperms)
  ci <- quantile(x = rdrb, probs = c(0.025,0.975))
  # fraction of randomized mean distances less than the observed mean
  f <- (sum( rdrb <= obsm)+1)/(nperms+1)
  out [i, ] <- c( nloc, obsm, nperms, ci, f)
}

out <- as.data.frame(out)
out$padj <- p.adjust(p = out$sig, method = 'holm')

out <- round(out, 4)
out$nperm <- round(out$nperm)

write.csv(x = out, file = fout, quote = F)



cat("Finished with a")
print(Sys.time()-tstart)
cat("\n")
