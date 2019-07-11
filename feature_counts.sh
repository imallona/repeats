#!/bin/bash
##
## part of the repeats scrnaseq profiler fast prototype
##
## Izaskun Mallona
## 11th July 2019

export SOFT=~/soft
export NTHREADS=20 ## number of threads to be used
export FEATURECOUNTS="$SOFT"/subread/subread-1.6.2-source/bin/featureCounts

export WD=/home/Shared_s3it/imallona/repeats_sc
export GTF="$WD"/data/rmsk.gtf
export PBMCS="$WD"/data/5k_pbmc_v3_possorted_genome_bam.bam

cd $WD

mkdir -p featurecounts

# duplicates are not ignored but multimappers are (and star reports plenty of them)
$FEATURECOUNTS \
    -T "$NTHREADS" \
    -t exon \
    -g gene_id \
    -a "$GTF" \
    -o  featurecounts/5k_pbmc_v3_possorted_genome_bam_rmsk.counts \
    "$PBMCS"
