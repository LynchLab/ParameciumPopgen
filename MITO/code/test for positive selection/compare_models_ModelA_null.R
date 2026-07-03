setwd('~/ASU Dropbox/Farhan ALI/Paramecium/intersp/final_315_samples/mito/dnds/')

fin <- 'modelA_null_lnL_values.csv'
fout <- 'modelA_null_lrt_out.csv'

## degree of freedoms
dof = 1
dof0 = 3
dof1 = 4

## significance threshold
alpha = 0.05

df <- read.csv( file = fin)

df <- df[ !apply( df, 1, anyNA), ]
n <- nrow(df)


get_pval <- function(l0,l1,dof){
  ld <- -2*(l0-l1)
  pchisq(q = ld, df = dof, lower.tail = F)
}

loglik_null <- as.numeric(df$null_lnL)
loglik_model <- as.numeric(df$ModelA_lnL)

df$p <- sapply( 1:n, function(i) get_pval(l0 = loglik_null[i], l1 = loglik_model[i], dof = dof))
df$adjp <- p.adjust( p = df$p, method = 'holm')

df$sig <- (df$adjp < alpha)

df <- df[ order(df$p, decreasing = F), ]

### AIC
df$aic0 <- 2*(dof0 - df$null_lnL)
df$aic1 <- 2*(dof1 - df$ModelA_lnL)

get_rell <- function(a1,a2){
  as <- c( a1, a2)
  mina <- min(as)
  bs <- exp( ( mina - as)/2)
  return(bs)
}

rels <- t( apply( df[ , c("aic0","aic1")], 1, function(x) get_rell(x[1],x[2])))
df$rel0 <- rels[,1]
write.csv(x = df, file = fout, quote = F, row.names = F)