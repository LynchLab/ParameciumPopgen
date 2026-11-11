fog <- 'single_burs_multicopy_aurelias_core_oma.tsv' # transcript id table of orthologs
fpm <- 'sp_prot_transcript_id_pfx.tsv'
indr <- 'data' # dir of species tables of genes' pnps
fout <- 'pnps_reldif_aure_burs_ncopies.csv'

dog <- read.csv(file = fog, sep = '\t', row.names = 1)
dpm <- read.csv(file = fpm, header = F, sep = '\t', row.names = 1)

## limit to species with nuclear genomes and population data
species <- c('primaurelia','biaurelia','tetraurelia','sexaurelia','bursaria')
nsp <- length(species)
# 4-letter abbreviations for species names
spshort <- sapply( strsplit(x = species, split = '', fixed = T), function(x) paste(x[1:4],collapse = ''))
# remove any extra bit in the names
colnames(dog) <- gsub("_.*", "", colnames(dog))
# burs is included as row names
dog <- dog[ , spshort[ -nsp]]

### convert ortholog table from wide to long format ####
genomes <- colnames(dog)
groups <- rownames(dog)

# list orthologs -> species -> vector of genes
# strings of genes are split on commas
og_sp_genes <- apply( dog, 1, function(x) strsplit( x, split = ',', fixed = T))
# attach group and species info to each vector of genes
lsog <- lapply( groups, function(x) lapply( 
  genomes, function(y) cbind( x, y, og_sp_genes[[x]][[y]])))
# convert the nested list into a long dataframe
ogl <- as.data.frame( do.call( rbind, unlist( lsog, recursive = F)))
rownames(ogl) <- NULL
colnames(ogl) <- c('group','species','protein')
## add transcript ids
tids <- dpm$V3
names(tids) <- dpm$V2
ppfx <- gsub('[0-9]+$', '', ogl$protein)
ogl$transcript <- sapply( 1:nrow(ogl), function(x)
  gsub( pattern = ppfx[x], replacement = tids[ppfx[x]], x = ogl$protein[x]))
########################################################

### get pN/pS of genes in each group from corresponding species ####
ogl$pnps <- NA

nsp <- length(genomes)

for(i in 1:nsp){
  # load species pnps data
  fin <- file.path( indr, paste0( species[i],".csv"))
  df <- read.csv( file = fin, row.names = 1, check.names = F)
  # fill in corresponding entries in the pooled data table
  sixs <- which( ogl$species == genomes[i])
  ogl[ sixs, 'pnps'] <- df[ ogl[ sixs, 'transcript'], 'pN/pS']
}

## pN/pS for corresponding bursaria homologs
burs_genes <- unique(ogl$group)
fin <- file.path( indr, 'bursaria.csv')
df <- read.csv(file = fin, row.names = 1, check.names = F)
burs_pnps <- df[ burs_genes, 'pN/pS']
names(burs_pnps) <- burs_genes
######################################################

### no. of homologs and mean pnps in each species for each group #####
ngps <- length(burs_genes)
sps <- unique(ogl$species)
nsp <- length(sps)
nmhom <- NULL

# for each ortholog group
for( i in 1:ngps){
  # subset the table
  sogl <- ogl[ ogl$group == burs_genes[i], ]
  # for each species
  for( j in 1:nsp){
    # further subset 
    ssogl <- sogl[ sogl$species == sps[j], ]
    # get no. of copies and average pnps over those copies
    entry <- c( burs_genes[i], sps[j], nrow(ssogl), mean(ssogl$pnps,na.rm = T))
    nmhom <- rbind( nmhom, entry)
  }
}
nmhom <- data.frame( group = nmhom[,1], species = nmhom[,2],
                     ngenes = as.integer(nmhom[,3]), 
                     avg_pnps = as.numeric(nmhom[,4]))

# attach pnps in bursaria
nmhom$burs_pnps <- burs_pnps[ nmhom$group]

# make sure that pnps are not NA or infinite
nmhom <- nmhom[ !is.na(nmhom$avg_pnps) & !is.na(nmhom$burs_pnps) &
                  !is.infinite(nmhom$avg_pnps) & !is.infinite(nmhom$burs_pnps), ]

#  bursaria pnps must be nonzero
nmhom <- nmhom[ nmhom$burs_pnps > 0,  ]

# calculate relative difference A - B / B
nmhom$difrelburs <- nmhom$avg_pnps / nmhom$burs_pnps - 1
# save 
write.csv(x = nmhom, file = fout, quote = F, row.names = F)
