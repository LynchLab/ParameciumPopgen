densityS_eq4_NY03 <- function(S,a,b){
  h = exp(S)-1
  t1 = b * exp(S) * S / h
  t2 = exp(-t1)
  t3 = h - S
  t4 = S * h * gamma(a)
  t1^a * t2 * t3 / t4
}

my_densityS <- function(S,a,b){
  h = exp(-S)-1
  w = -S/h
  dwdS <- -(h+S*exp(-S))/h^2
  fw <- dgamma(x = w,shape = a,rate = b)
  return(fw*dwdS)
}

cdf_S <- function(S,a,b){
  h = exp(-S)-1
  w = -S/h
  pgamma(q = w,shape = a,rate = b,lower.tail = F)
}

quantile_S <- function(a,b,p=c(0.025,0.975)){
  q <- qgamma(p = p, shape = a,rate = b)
  h = exp(-q)-1
  w = -q/h
  w
}

exp_integrand_plus <- function(x,a,b){
  ifelse(test = x == 0, yes = 0, no = x*my_densityS(x,a,b))
}


## log sum exp
log_sum_exp <- function(logv) {
  m <- max(logv, na.rm = TRUE)
  if (!is.finite(m)) return(-Inf)
  return(m + log(sum(exp(logv - m), na.rm = TRUE)))
}

exp_integrand_neg <- function(t,a,b){
  # replacing S w -t
  h = exp(t)-1
  # log h(-t)
  q = t + log(1-exp(-t))
  em <- exp( -b*t*exp(-q) + a*log(b*t) + t - q*(1+a))*(1-t-exp(-t))
  ifelse(test = t == 0, yes = 0, no = em)
}

exp_S <- function(a,b){
  eplus <- integrate( f = exp_integrand_plus, lower = 0, upper = Inf, a = a, b = b)$value
  emins <- integrate( f = exp_integrand_neg, lower = 0, upper = Inf, a = a, b = b)$value
  out <- eplus + emins/gamma(a)
  return(out)
}


### examples ######
a = .373
b = .523

Svec <- seq(-20,20,by=0.1)

dS <- densityS_eq4_NY03(S = Svec,a = a,b = b)

mydS <- my_densityS(S = Svec,a = a,b = b)

plot( Svec, mydS, type = 'l')
abline( v = c(-1,0,1),lty=3)

eS <- exp_S(a,b)
abline( v = eS)

############################
setwd('~/ASU Dropbox/Farhan ALI/Paramecium/intersp/final_315_samples/mito/dnds')
fin <- 'M5/M5_gamma_params.csv'
fout <- 'M5/gene_S_drbs.pdf'
df <- read.csv(file = fin, row.names = 1)
ng <- nrow(df)

Svec <- seq(-20,20,0.01)

# initialize vector of gene-wise expectation of S
gexps <- rep( NA, ng)
names(gexps) <- rownames(df)

pdf(fout)
par(mfrow = c(4,4), mar = c(3,3,2,2))

for(i in 1:ng){
  gdens <- my_densityS(S = Svec, a = df[i,'a'], b = df[i,'b'])
  # total dens
  td <- sum(gdens,na.rm=T)
  # mean S
  mS <- round( exp_S(a = df[i,'a'],b = df[i,'b']),1)
  gexps[i] <- mS
  # proportion of sites w S > 0
  p <- round(cdf_S(S = 1,a = df[i,'a'], b = df[i,'b']),3)
  # 95 % quantiles
  ql <- round(cdf_S(S = 1,a = df[i,'a'], b = df[i,'b']),3)
  plot( Svec, gdens, type = 'l', xlab = '', ylab = '')
  abline( v = c(-1,0,1),lty=3,col='grey40')
  abline( v = mS, lty = 2)
  mtext(text = rownames(df)[i], adj = 0.9, line = -1, cex = 0.7)
  mtext(text = bquote( paste( bar(italic(S)), ' = ', .(mS))), side = 3, adj = 0, cex = 0.5)
  mtext(text = bquote( paste( italic(p)[S>1], ' = ', .(p))), side = 3, adj = 1, cex = 0.7)
}

dev.off()

### correlation with overall omega ###
fom <- 'overall/overall_dnds_gene.csv'
dom <- read.csv( file = fom, row.names = 1)

genes <- intersect( rownames(dom), names(gexps))
gomms <- data.frame( om = dom[genes,], es = gexps[genes])

plot( es ~ om, data = gomms)

### correlation with pairwise omega #######
fpm <- 'pairwise/'

## which one correlates better with pN/pS?
fpnps <- '../pnps/pnps_mean_all_species.csv'
dpnps <- read.csv( file = fpnps)
# filter out low number of syn differences
dpnps <- dpnps[ dpnps$Mean_Ds > 0, ]
# group by genes
gene_pnps <- split( dpnps$pN_over_pS, dpnps$Gene)
gene_pnps_summ <- as.data.frame( t( sapply( gene_pnps, summary)))
# add number of points for a gene
gene_pnps_summ$n <- sapply( gene_pnps, length)
# reorder by median
gene_pnps_summ <- gene_pnps_summ[ order(gene_pnps_summ$Median), ]
gene_pnps_summ

#### atp 9 #################
f1 <- 'overall/atp9_branch_length_estimates.txt'
f2 <- 'M5/atp9_branch_length_estimates.txt'

df1 <- read.csv(file = f1,sep = ',',check.names = F, row.names = 1)
df2 <- read.csv(file = f2,sep = ',',check.names = F, row.names = 1)

head(df1)
head(df2)

bnames <- intersect( rownames(df1), rownames(df2))
df1 <- df1[ bnames, ]
df2 <- df2[ bnames, ]

plot( df1$t, df2$t, log = 'xy')
abline(0,1)

#### ymf 64 #################
f1 <- 'overall/ymf64_bls.txt'
f2 <- 'M5/ymf64_bls.txt'

df1 <- read.csv(file = f1,sep = ',',check.names = F, row.names = 1)
df2 <- read.csv(file = f2,sep = ',',check.names = F, row.names = 1)

head(df1)
head(df2)

bnames <- intersect( rownames(df1), rownames(df2))
df1 <- df1[ bnames, ]
df2 <- df2[ bnames, ]

plot( df1$t, df2$t)
abline(0,1)

### numerical methods to convert omega to S ###########
rootfunc <- function(x,w){
  x/(1-exp(-x))-w
}

gx <- rootfunc(Svec,0.1)
plot( Svec, gx, type = 'l', xlim = c(-10,5), ylim = c(-1,5))
abline(h = 0, lty = 3)

nrest_ob <- uniroot(f = rootfunc, interval = c(-20,20), w = 0.1)
estS <- nrest_ob$root
abline(v = estS, lty = 2)

#### with beta distribution on omega ####################
beta_densityS <- function(S,p,q){
  h = exp(-S)-1
  w = -S/h
  dwdS <- -(h+S*exp(-S))/h^2
  fw <- dbeta(x = w, shape1 = p, shape2 = q)
  return(fw*dwdS)
}

beta_cdfS <- function(S,p,q){
  h = exp(-S)-1
  w = -S/h
  pbeta( q = w, shape1 = p, shape2 = q, lower.tail = F)
}


fin <- 'M8/M8_gene_params.csv'
fout <- 'M8/gene_S_beta_density.pdf'
df <- read.csv(file = fin, row.names = 1)
df <- df[ df$p0 > 0.999, ]
ng <- nrow(df)

Svec <- seq(-20,5,0.01)
neg_Svec <- Svec[ Svec < 0]
n_neg_svec <- length(neg_Svec)

# initialize vector of gene-wise expectation of S
gexps <- rep( NA, ng)
names(gexps) <- rownames(df)

pdf(fout)
par(mfrow = c(2,2), mar = c(3,3,2,2))

for(i in 1:ng){
  gdens <- beta_densityS(S = Svec, p = df[i,'p'], q = df[i,'q'])
  neg_gdens <- gdens[1:n_neg_svec]
  modeS <- neg_Svec[ which( neg_gdens == max(neg_gdens))]
  # proportion of sites w S > 1
  p <- round( beta_cdfS(S = 1, p = df[i,'p'], q = df[i,'q']),3)
  plot( Svec, gdens, type = 'l', xlab = '', ylab = '', xaxt = 'n')
  abline( v = c(-1,0,1),lty=3,col='grey40')
  axis(side = 1, at = seq(-20,5,by=5))
  abline( v = modeS, lty = 2)
  mtext(text = rownames(df)[i], adj = 0.9, line = -1, cex = 0.7)
  mtext(text = bquote( paste( bar(italic(S)), ' = ', .(modeS))), side = 3, adj = 0, cex = 0.5)
  #mtext(text = bquote( paste( italic(p)[S>1], ' = ', .(p))), side = 3, adj = 1, cex = 0.7)
}

dev.off()
