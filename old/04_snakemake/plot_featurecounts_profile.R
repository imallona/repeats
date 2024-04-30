#!/usr/bin/env R

## gets the featurecounts stats (input) and generate a barchart of annotated/unannotated etc
##
## if  multiple cells present, plots four randomly sampled
##
## 21 May 2020

## install.packages('Cairo', repos ="https://cloud.r-project.org")

suppressPackageStartupMessages({
    library("optparse")
    library('data.table')
    library('reshape2')
    library('ggplot2')
    library('Cairo')
})


option_list = list(
    make_option(c("-s", '--summary'),
                type = 'character', default = NULL,
                help = 'summary file as produced by featurecounts'),
    make_option(c("-o", '--output'),
                type = 'character', default = NULL,
                help = 'png output file'),
    make_option(c("-i", '--identifier'),
                type = 'character', default = NULL,
                help = 'gtf identifier'))
    
options(bitmapType='cairo')

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

## opt$summary <- '/home/imallona/repeats_sc/runs/ding_celseq2/bowtie_repeatome/ding_celseq2SRR9169177_bowtie_repeats.counts.summary'
## opt$id <- 'test'
## opt$summary <- '/home/imallona/repeats_sc/runs/pbmc_10k_v3/count_repeats_on_cellranger_standard/pbmc_10k_v3_repeats.counts.summary'
## opt$output <- '/home/imallona/repeats_sc/runs/pbmc_10k_v3/count_repeats_on_cellranger_standard/transcriptome_repeats.counts.summary.png'
## opt$identifier <- 'cellranger_repeats'



fd <- data.table::fread(opt$summary, header = TRUE)

colnames(fd) <- basename(colnames(fd))
if (ncol(fd) > 2) {
    set.seed(2574)
    idx <- sample(x = 2:ncol(fd), size = 4, replace = FALSE)
    fd <- as.data.frame(fd)[,c(1, idx)]
}

pct <- fd

pct <- as.data.frame(apply(pct[,2:ncol(pct)], 2, prop.table))
pct$Status <- fd$Status

pct <- melt(pct, id.vars=c("Status"))
fd <- melt(fd, id.vars=c("Status"))


pct$counts <- fd$value
colnames(pct) <- c('Status', 'bam', 'pct', 'count')
pct$pct <- pct$pct * 100

p <- ggplot(pct, aes(x=Status, y = pct)) +
    geom_bar(stat="identity") +
    facet_wrap(~bam) +
    geom_text(aes(label=sprintf('%.2g', count)), size=3,
              vjust = -0.5) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ylab ('read assignments (%)') +
    ggtitle(sprintf('GTF %s', opt$id))

## ggsave(p, filename = sprintf('%s_barplot.png', opt$output), device = 'png')

write.csv(pct, file = sprintf('%s_barplot.csv', opt$output))

ggsave(p, filename = opt$output, device = 'png')
