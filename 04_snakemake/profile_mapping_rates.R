#!/usr/bin/env R

## gets the featurecounts stats (input) and generate a barchart of annotated/unannotated etc
## as well as the alevin's and cellranger
##
## 17 Feb 2021


suppressPackageStartupMessages({
    library("optparse")
    library('data.table')
    library('reshape2')
    library('ggplot2')
    library('Cairo')
    library('jsonlite')
})

option_list = list(
    make_option(c("-p", '--path'),
                type = 'character', default = NULL,
                help = 'path with logs'),
    make_option(c("-o", '--output'),
                type = 'character', default = NULL,
                help = 'final ouput file (full path), csv'),
    make_option(c("-i", '--identifier'),
                type = 'character', default = NULL,
                help = 'run identifier'))
    

process_featurecounts_logs <- function(fns, pattern = NA) {
    if (!is.na(pattern))
        fns <- grep(pattern, fns, value = TRUE)
    d <- list()

    for (fn in fns) {
        
        tmp <- rowSums(data.frame(data.table::fread(fn, header = TRUE), row.names = 1))
        
        ## colnames(fd) <- basename(colnames(fd))

        fd <- data.frame(assigned = tmp['Assigned'],
                         total = tmp['Unassigned_Unmapped'] +tmp['Unassigned_Read_Type'] +
                             tmp['Unassigned_Singleton'] +tmp['Unassigned_MappingQuality'] +
                             tmp['Unassigned_Chimera'] +tmp['Unassigned_FragmentLength'] +
                             tmp['Unassigned_Duplicate'] +tmp['Unassigned_MultiMapping'] +
                             tmp['Unassigned_Secondary'] +tmp['Unassigned_NonSplit'] +
                             tmp['Unassigned_NoFeatures'] +tmp['Unassigned_Overlapping_Length'] +
                             tmp['Unassigned_Ambiguity'])

        d[[fn]] <- fd
    }
    ## d
    prop.table(Reduce(`+`, d))    
}

## fns <- fns$alevins
process_alevin_logs <- function(fns, pattern = NA) {
    if (!is.na(pattern))
        fns <- grep(pattern, fns, value = TRUE)
    d <- list()

    ## the mapping rate is just $num_mapped / $num_processed
    for (fn in fns) {
        tmp <- read_json(fn)
        d[[fn]] <- setNames(c(tmp$num_mapped, tmp$num_processed), c('mapped', 'total'))
        ## d[[fn]] <- setNames(read_json(fn)$mapping_rate, 'mapping_rate')
    }
    prop.table(Reduce(`+`, d))  
}

## fns <- fns$cellrangers
process_cellranger_logs <- function(fns, pattern = NA) {
    if (!is.na(pattern))
        fns <- grep(pattern, fns, value = TRUE)
    d <- list()
    for (fn in fns) {
        fd <- data.table::fread(fn, header = TRUE)
        d[[fn]] <- setNames(c(
            as.numeric(gsub('%', '' ,fd$`Reads Mapped Confidently to Transcriptome`))/100 * as.numeric(gsub(',', '', fd$`Number of Reads`)),
            as.numeric(gsub(',', '', fd$`Number of Reads`))),
                            c('mapped', 'total'))
    }
    prop.table(Reduce(`+`, d))  
}

process_logs <- function(fns, software, file_pattern = NA) {
    if (software == 'cellranger') {
        process_cellranger_logs(fns = fns, pattern = file_pattern)
    }
    else if (software == 'featurecounts') {
        process_featurecounts_logs(fns = fns, pattern = file_pattern)
    }
    else if (software == 'alevin') {
        process_alevin_logs(fns = fns, pattern = file_pattern)
    }
    else {
        stop('software not recognised')
    }
}

options(bitmapType='cairo')

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (FALSE) {
    opt <- list()
    opt$path <- '/home/imallona/repeats_sc/runs/pbmc8k'
    opt$identifier <- 'test'
    opt$output <- '~/tmp/test.csv'
}

## opt$identifier
fns <- list()
fns$featurecounts <- file.path(
    opt$identifier,
    list.files(opt$identifier, pattern = "*counts\\.summary$?", recursive = TRUE))

fns$alevin <- file.path(
    opt$identifier,
    list.files(opt$identifier, pattern = "^meta_info.json", recursive = TRUE))

fns$cellranger <-  file.path(
    opt$identifier,
    list.files(opt$identifier, pattern = "^metrics_summary.csv$", recursive = TRUE))


## fns$alevin
## grep('*alevin.*repeats*', fns$alevin, value =TRUE
     )

res <- data.frame(pattern = c('cellranger_repeats',
                              'cellranger_standard',
                              'alevin.*repeats*',
                              'alevin.*genes*',
                              'bowtie_repeatome/multimappers',
                              'bowtie_repeatome/unique_reads',
                              'count_repeats_on_cellranger_standard_not_overlapping_genes/multimappers',
                              'count_repeats_on_cellranger_standard_not_overlapping_genes/unique_reads',
                              'count_repeats_on_cellranger_standard/multimappers',
                              'count_repeats_on_cellranger_standard/unique_reads'))

    
res$dataset <- opt$identifier
res$analyte <- c(rep(c('repeats', 'genes'), 2), rep('repeats', 6)) 
res$mapper  <- c(rep('cellranger', 2), rep('alevin', 2), rep('bowtie', 2), rep('cellranger', 4))
res$counter <- c(rep('cellranger', 2), rep('alevin', 2), rep('featurecounts', 6))
res$overlap <- c('repeats', 'genes', 'decoy', 'genes', 'repeats', 'repeats', 'repeats_only', 'repeats_only', 'repeats', 'repeats' )
res$multimapping <- c(rep('cellranger', 2), rep('alevin', 2), rep(c('multi', 'unique'),3))
res$mapped <- NA
res$total <- NA

for (i in 1:nrow(res)) {
    software <- res[i, 'counter']
    
    pattern <- res[i, 'pattern']
    tryCatch({
        res[i, c('mapped', 'total')] <- process_logs(fns = fns[[software]], software = software,
                                                     file_pattern = pattern)
        }, error = function(x) print(x))
    
}


write.csv(pct, file =  opt$output)
