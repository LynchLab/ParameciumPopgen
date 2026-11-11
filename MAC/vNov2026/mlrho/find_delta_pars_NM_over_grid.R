set.seed(112211)

args <- commandArgs(trailingOnly = T)
indr <- args[1]
fout <- args[2]
fsfx <- args[3]

mindists <- args[4] # minimum distance threshold

tstart <- Sys.time()

# bounds on each parameter on the original scale
 lims <- data.frame( theta = c( 1e-04, 0.5), 
                     rho1 = c( 1e-04, 0.5),
                     x = c( 0.001, 0.999), 
                     L = c( 10, 1e04))


######## FUNCTIONS #################################################################
## function to filter grid based on objective
filter_grid <- function( grid, data, objective, maxer=1){
  ## check objective's value on each grid point
  grid$er <- apply( grid, 1, function(x)
    objective( pars = x, data = data))
  
  # find the best - minimum error
  miner <- min(grid$er)
  bestx <- which( grid$er == miner)
  
  # error relative to the best
  grid$reler <- grid$er / miner - 1
  
  # keep grid points within relative maxer
  grid <- grid[ grid$reler <= maxer, , drop = F]
  return(grid)
}

## spacings on each axis of a uniform grid
get_spacing <- function(grid){
  spacings <- apply( grid, 2, function(x) diff( sort( unique(x))[1:2]))
  # for axis reduced to a point, set h to 1 
  # so that the axis remains the same on the normalized grid
  spacings[ is.na(spacings)] = 1
  return(spacings)
}

## find first K best grid points that are at least a given width away
thin_grid <- function( grid, hvec, bestK=5, minmd=2){
  # no. of dimensions
  ndim = length(hvec)
  # no. of grid points
  ngp <- nrow(grid)
  if(ngp <= 1 || bestK == 1)
    return(grid)
  
  # order grid, best first
  grid <- grid[ order(grid$reler), ]
  
  # normalize grid by dividing each axis by its spacing
  normgrid <- t(t( grid[,1:ndim])/hvec)
  
  # initialize selected indices with the first point
  selxs <- 1
  
  # pick k best
  k = 1
  while(k < bestK){
    # last included index
    lx <- selxs[k]
    # check that last index is not at the last grid point
    if( lx == ngp) break
    
    # remaining indices to check
    cixs <- (lx+1):ngp
    
    # which ones satisfy gap threshold to all collected values
    xpass <- which( sapply( cixs, function(j)  
      min( apply( normgrid[ selxs, , drop = F], 1, function(x)
        sum( abs( x - normgrid[j,])))) >= minmd))
    
    # check if any does
    if( length(xpass) == 0) break
    
    # if some do, then add the first to selected indices
    selxs <- c( selxs, cixs[ xpass[1]])
    # update k
    k <- length(selxs)
  }
  
  grid[ selxs, , drop = F]
}

## main equation for expected delta as a function of distance and parameters
expected_delta <- function(d, theta, rho1, x, L) {
  # Eq. 13a: distance-dependent recombination rate (population-scaled)
  rho_d <- rho1 * (x * d + (1 - x) * L * (1 - exp(-d / L)))
  # Eq. 5d: A
  A <- 9 + 6.5*rho_d + 0.5*rho_d^2 + 12*theta*(4+5*theta+2*theta^2) +
    rho_d*theta*(17 + rho_d + 8*theta)
  # Eq. 6
  theta*(1+2*theta)*(18+rho_d)/(2*A*(1+theta))
}

## Objective function to minimize: Mean squared deviation in log space
objective <- function(pars, data) {
  # fed-parameters are log-transformed
  # so back-transformation ensures positivity
  theta <- exp(pars[1])  
  rho1  <- exp(pars[2])
  x     <- plogis(pars[3])  # logistic to enforce x in (0,1)
  L     <- exp(pars[4])
  # corresponding delta
  delta_pred <- expected_delta( d = data[,1], theta, rho1, x, L)
  
  # guard against non-positive predictions
  if (any(delta_pred <= 0)) return(1e10)
  # average squared difference in log space
  mean((log( data[,2]) - log(delta_pred))^2)
}

## Nelder-Mead optimization on a single grid point
fit_lsq <- function( init, data){
 
  fit <- optim(
    par     = init,
    fn      = objective,
    data    = data,
    method  = "Nelder-Mead",
    control = list(maxit = 10000, reltol = 1e-10)
  )
  return(fit)
}


## function to generate a coarse grid and filter given data
generate_grid <- function( lims, data, fun, bestK, nstep=5,maxer=1){
  ####################
  # lims - dataframe with 2 rows and m columns for bounds on m parameters to estimate
  # data - dataframe or matrix with 2 columns - (X,Y)
  # fun - objective function to minimize the value of over parameters. eg. least squares
  # bestK - K initial best guesses from the coarse grid to do a local search over
  # nstep - number of steps on either side of the center point
  # nuz - number of units in the initial grid to zoom over
  ####################
  # no. of dimensions
  ndim <- ncol(lims)
  ## generate initial coarse grid
  parvals <- lapply( 1:ndim, function(i) 
    seq( lims[1,i], lims[2,i], length.out = 2*nstep+1))
  
  grid <- expand.grid( parvals)
  cat(sprintf("size of initial coarse grid = %i x %i\n",nrow(grid),ncol(grid)))
  
  ## filter grid
  fgrid <- filter_grid( grid[ , 1:ndim], data = data, objective = fun, maxer = maxer)
  cat(sprintf("grid points left post-filtering by error = %i\n",nrow(fgrid)))
  
  ## spacings on each axis on the original grid
  ogspacings <- get_spacing( grid)
  
  ## thin out grid to K best points
  tgrid <- thin_grid(grid = fgrid, hvec = ogspacings, bestK = bestK)
  return(tgrid)
}

#### main algorithm #########################################
find_delta_pars <- function( obs_data, gridlims, mind=mindists){
  print(Sys.time()-tstart)
  
  # filter distances
  obs_data <- obs_data[ obs_data$dist >= mind, ]
  # initial generation of grid + filtering based on data
  cands <- generate_grid( gridlims, obs_data, objective, bestK = 5)
  
  # NM run on each
  tout <- apply( cands[ , 1:4, drop = F], 1, function(x) fit_lsq(init = x, data = obs_data))
  tout <-  t( sapply( tout, function(x) c(x$par,x$value)))
  
  # check that all estimated parameters are within limits
  ilxs <- apply( tout[,1:4,drop=F], 1, function(x) all( 
    x > gridlims[1,] & x < gridlims[2,]))
  
  if( sum(ilxs) == 0){
    cat("\n\tEstimated values exceed parameter limits!\n")
    badout <- tout[ !ilxs, ]
    badout[,1] <- exp(badout[,1])
    badout[,2] <- exp(badout[,2])
    badout[,4] <- exp(badout[,4])
    badout[,3] <- plogis(badout[,3])
    print(badout)
    cat("\n")
    return(NULL)
  }
  
  tout <- tout[ ilxs, , drop = F]
  
  # best set of parameters among them
  best <- tout[ tout[,5] == min(tout[,5]), , drop = F]
  best <- as.numeric(best[ 1, ])
  # back-transform parameters
  bbest <- c( exp(best[1]), exp(best[2]), plogis(best[3]), exp(best[4]))
  
  # output
  pout <- c( bbest, best[5])
  return( pout)
}

##############################################################################

### generate search grid for initial values  ######

## bounds on transformed scale
limstr <- data.frame( theta = log(lims$theta), rho1 = log(lims$rho1),
                      x = qlogis(lims$x), L = log(lims$L))

#################################################################


### DATA ####################################################
fls <- list.files(path = indr, pattern = fsfx)
nf <- length(fls)

# list of dataframes
lsdfs <- lapply( fls, function(f) 
  parse_delta_dist( file.path( indr, f))[,c("dist","delta")])

# sample distances past 1000 in intervals of 10
subdists <- c( mindists:999, seq(1000,5000,10))
## keep distances in multiples of 3
subdists <- subdists[ subdists%%3 == 0]

lsdfs <- lapply( lsdfs, function(X) X[ X$dist %in% subdists, ])

# initialize output table
out <- as.data.frame( matrix( NA, 0, 6, dimnames = list(
  NULL, c( 'isolate', 'theta', 'rho1', 'x', 'Lbar', 'error'))))
write.csv( x = out, file = fout, quote = F, row.names = F)

# run parameter search for each dataset
for( i in 1:nf){
  cat(sprintf("\ndataset # = %i, %s\n",i,fls[i]))
  fpfx <- gsub( fsfx, "", fls[i], fixed = T)
  
  tout <- tryCatch( find_delta_pars(obs_data = lsdfs[[i]], gridlims = limstr),
            error = function(e) {message("Error"); return(NULL)})
  if(is.null(tout)) next
  cat("\n")
  linedf <- data.frame( isolate = fpfx, theta = tout[1], rho1 = tout[2], 
                        x = tout[3], Lbar = tout[4], error = tout[5])
  write.table( x = linedf, file = fout, sep = ',', append = T, 
               quote = F, col.names = F, row.names = F)
}
