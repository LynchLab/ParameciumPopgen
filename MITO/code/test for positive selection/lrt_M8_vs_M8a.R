setwd('~/ASU Dropbox/Farhan ALI/Paramecium/intersp/final_315_samples/mito/dnds/')

fin <- 'm8_vs_m8a_lrt_raw.tsv'

dof0 = 3
dof1 = 4
dof = dof1 - dof0

alpha = 0.01

df <- read.csv( file = fin, row.names = 1, sep = '\t')[,1:2]

get_pval <- function(l0,l1,dof){
  ld <- -2*(l0-l1)
  pchisq(q = ld, df = dof, lower.tail = F)
}

colnames(df) <- c('M8a','M8')

df$p <- apply( df, 1, function(x) get_pval(l0 = x[1], l1 = x[2], dof = dof))
df$adjp <- p.adjust( p = df$p, method = 'holm')

df$sig <- (df$adjp < alpha)

df <- df[ order(df$adjp), ]

### AIC
df$aic0 <- 2*(dof0 - df$M8a)
df$aic1 <- 2*(dof1 - df$M8)

get_rell <- function(a1,a2){
  as <- c( a1, a2)
  mina <- min(as)
  bs <- exp( ( mina - as)/2)
  return(bs)
}

rels <- t( apply( df[ , c("aic0","aic1")], 1, function(x) get_rell(x[1],x[2])))
df$rel0 <- rels[,1]