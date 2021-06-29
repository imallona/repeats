#!/usr/bin/env R
##
## Parses a rmsk mysql query to generate a gtf
##
## Writes stuff to stdout (line by line)
##
## 14 June 2021
## 

suppressPackageStartupMessages({
    library("optparse")
})


option_list = list(
    make_option(c("-r", '--rmsk'),
                type = 'character', default = NULL,
                help = 'rmsk csv dump'),
    make_option(c("-s", '--source'),
                type = 'character', default = NULL,
                help = 'source field for column 2 in gtf'))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

## print(sprintf('processing GTF %s', opt$gtf))

## opt <- list(rmsk = 'tmp.dump.gz')
con <- file(opt$rmsk, open = "r")


while (TRUE) {
    line <- readLines(con, n = 1)
    if (length(line) == 0) {
        break
    }

    line <- strsplit(line, split = '\t')[[1]]
    seqname <- line[1]
    source <- sprintf('%s_rmsk', opt$source)
    feature <- 'exon'
    start <- line[2]
    end <- line[3]
    score <- line[4]
    strand <- line[5]
    frame <- "0"
    transcript_id <- sprintf('%s_%s_%s_%s_%s', line[6], seqname, start, end, strand)
    
    name_id <- sprintf('gene_id "%s"; ', line[6])
    family_id <- sprintf('family_id "%s"; ', line[8])
    
    class_id <- sprintf('class_id "%s"; ', line[7])

    attribute <- paste(c(name_id, sprintf('transcript_id "%s"; ', transcript_id),
                         family_id, class_id), collapse = '')

    cat(sprintf('%s\n', paste(c(seqname, source, feature,
                start, end,
                score, strand, frame, attribute),
                collapse = '\t')))
}

close(con)
