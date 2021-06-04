#!/usr/bin/env R
##
## Shao's TE assembly simulation
##

library(data.table)

saf <- '~/repeats_sc/annotation/shao/shao_GRCh38_te.saf'
truth <- '~/repeats_sc/data/sim_5k_pbmc_v3/truth.tsv.gz'
counts <- '~/repeats_sc/runs/sim_5k_pbmc_v3/zumis/sim_5k_zumis_count_dt.txt'

d <- list()
d$saf <- fread(saf)
d$truth <- fread(truth)
d$counts <- fread(counts)

sapply(d, colnames)

sapply(d, nrow)

merged <- merge(d$counts, d$saf, by.x= 'feature', by.y = 'GeneID', all = TRUE)
dim(merged)
head(merged)

head(d$truth)

d$truth$found <- NA

d$truth$chr <- gsub(".*::(.+):([0-9]+)-([0-9]+)", "\\1", d$truth$locus)
d$truth$start <- gsub("(.+):([0-9]+)-([0-9]+)", "\\2", d$truth$locus)
d$truth$end <- gsub("(.+):([0-9]+)-([0-9]+)", "\\3", d$truth$locus)

head(d$truth)

## merged <- as.data.frame(merged)
## for (i in 1:nrow(d$counts)) {
##     feature <- d$counts$feature[i]
##     tmp <- merged[merged$feature == feature,]
    
    
    
## }
