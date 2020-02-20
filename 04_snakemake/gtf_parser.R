#!/usr/bin/env R
##
## Parses a gz-compressed GTF file to replace the third column by the
##  'transcript' field from column 9. Assumes all records are exons.
## This happens in rmsk-based GTFs.
##
## Writes stuff to stdout (line by line)
##
## 20 Feb 2020
## 

suppressPackageStartupMessages({
    library("optparse")
})


option_list = list(
    make_option(c("-g", '--gtf'),
                type = 'character', default = NULL,
                help = 'gz-compressed gtf to parse'))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

## print(sprintf('processing GTF %s', opt$gtf))

con <- gzfile(opt$gtf, "r")

while (TRUE) {
    line <- readLines(con, n = 1)
    if (length(line) == 0) {
        break
    }

    tmp <- strsplit(line, split = '\t')[[1]]

    tmp[3] <- sub('.*transcript_id\ "(.*)\"; family_id.*', "\\1", tmp[[9]])
    tmp[9] <- gsub('\"', '"', tmp[9])
    
    cat(sprintf('%s\n', paste(tmp, collapse = '\t')))
}

close(con)

