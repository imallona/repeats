#!/usr/bin/env R

# library(gtools)
library(DNABarcodes)

options <- commandArgs(trailingOnly = TRUE)
options



## LENGTH <- 92
## NTHREADS <- 20
## QUAL_SYMBOL <- 'F'
## # get rmsk coordinates
## GTF <- '~/repeats_sc/annotation/mm10_rmsk_TE.gtf.gz'

## GENOME_FASTA <- '~/repeats_sc/annotation/Mus_musculus.GRCm38.dna.primary_assembly.fa'

## BEDTOOLS <- '~/soft/bedtools/bedtools-2.29.2/bin/bedtools'


## I1_ORIGIN <- '~/repeats_sc/data/5k_pbmc_v3/5k_pbmc_v3_fastqs/5k_pbmc_v3_S1_L002_I1_001.fastq.gz'
## R1_ORIGIN <- '~/repeats_sc/data/5k_pbmc_v3/5k_pbmc_v3_fastqs/5k_pbmc_v3_S1_L002_R1_001.fastq.gz'
## R2_SIM_ORIGIN <- '~/repeats_sc/data/sim/5k_pbmc_v3_sim_S1_L002_R2_001.fastq'


## I1_SIM <- '~/repeats_sc/data/sim/5k_pbmc_v3_sim_S1_L002_I1_001.fastq'
## R1_SIM <- '~/repeats_sc/data/sim/5k_pbmc_v3_sim_S1_L002_R1_001.fastq'

R1_ORIGIN <- options[1]
R2_SIM_ORIGIN <- options[2]
I1_SIM <- options[3]
R1_SIM <- options[4]
NTHREADS <- as.numeric(options[5])
CB <- as.numeric(options[6])
UMI <- as.numeric(options[7])
NUM_CELLS <- as.numeric(options[8])
TRUTH_COUNT_TABLE <- options[9]

## CB <- 16
## UMI <- 12

## NUM_CELLS <- 1000

## simulate R2 file, one repeat per line, 92 nt-long lines

## cmd <- sprintf('simulate_cdna %s %s %s %s %s %s %s',
##                LENGTH,
##                QUAL_SYMBOL,
##                GTF,
##                GENOME_FASTA,
##                NTHREADS,
##                BEDTOOLS, R2_SIM_ORIGIN)

## print(cmd) # run this!

## to read the real R1s, so we can get some realistic CBs
true1r <- file(R1_ORIGIN, "r")

tmp <- readLines(true1r, n = 1e6)
tmp <- tmp[seq (2, length(tmp),4)]
stopifnot(all(nchar(tmp) == 28))

cb_freq <- as.data.frame(table(substr(x = tmp, start = 1, stop = CB)))
cbs <- as.character(head(cb_freq[order(cb_freq$Freq, decreasing = TRUE), 'Var1'],
           NUM_CELLS))

close(true1r)

## UMIs are just simulated, one each

umis <- create.dnabarcodes(n = UMI,
                           dist = 4,
                           heuristic = 'conway',
                           ## heuristic = "ashlock",
                           cores = NTHREADS)


## to read the simulated R2s, including the read names
sim2r <- file(R2_SIM_ORIGIN, "r")

## to write the simulated ones (I1)
sim1iw <- file(I1_SIM, "w")
sim1rw <- file(R1_SIM, "w")
truthw <- file(TRUTH_COUNT_TABLE, "w")

template_i1 <- '%s\nCATTAGCG\n+\nFFFFFFFF'
template_r1 <- '%s\n%s%s\n+\nFFFFFFFFFFFFFFFFFFFFFFFFFFFF'

writeLines(sprintf('%s\t%s\t%s\t%s', 'cell', 'umi', 'locus', 'count'), con = truthw)

curr_cell <- 1
curr_umi <- 1

while (TRUE) {    
    stanza <- readLines(sim2r, n = 4)
    if ( length(stanza) == 0 ) {
      break
    }
    id <- stanza[1]

    ## pick a CB randomly
    curr_cell <- sample(1:NUM_CELLS, 1)
    
    writeLines(sprintf(template_i1, id), con = sim1iw)
    writeLines(sprintf(template_r1, id,
                       cbs[curr_cell],
                       umis[curr_umi]), con = sim1rw)

    #@ we only simulate 1 read per cell, that's the last column
    writeLines(sprintf('%s\t%s\t%s\t%s', cbs[curr_cell], umis[curr_umi],
                       gsub('@', '', id), 1), con = truthw)

    ## get an UMI from the list
    curr_umi <- curr_umi +1
    
    ## if we exceed the number of umis, they're recycled (they'll likely map to different places)
    if (curr_umi > length(umis)) {
        curr_umi <- 1
    }
    
    if (curr_umi %% 1e4 == 0) cat('.') 
}

close(sim1iw)
close(sim1rw)
close(sim2r)
close(truthw)
