#!/usr/bin/env R
##
## Reads a fastq and mutates each record with as many mutations as desired
##  that are inserted randomly

options <- commandArgs(trailingOnly = TRUE)
options

FASTQ_INPUT <- options[1]
NUM_VARIANTS <- options[2]
FASTQ_OUTPUT<- options[3]

alphabet <- c('A', 'C', 'T', 'G')
add_mutation <- function(x, idx, alphabet){
    splitted <- strsplit(x, '')[[1]]
    for (i in idx) {
        splitted[i] <- sample(x = setdiff(alphabet, splitted[i]), size = 1)
    }    
    return(paste(splitted, collapse = ''))
}


origin <- file(FASTQ_INPUT, "r")
out <- file(FASTQ_OUTPUT, "w")

while (TRUE) {    
    stanza <- readLines(origin, n = 4)
    if ( length(stanza) == 0 ) {
      break
    }
    cdna <- stanza[2]

    pos <- sample(x = 1:length(cdna), size = NUM_VARIANTS, replace = FALSE)
    cdna <- add_mutation(x = cdna, idx = pos, alphabet = alphabet)

    writeLines(sprintf('%s\t%s\t%s\t%s', stanza[1], cdna, stanza[3], stanza[4]),
               out)
}

close(origin)
close(out)
