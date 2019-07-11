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

# # duplicates are not ignored but multimappers are (and star reports plenty of them)
# $FEATURECOUNTS \
#     -T "$NTHREADS" \
#     -t exon \
#     -g gene_id \
#     --byReadGroup \
#     -a "$GTF" \
#     -o  featurecounts/5k_pbmc_v3_possorted_genome_bam_rmsk_rg.counts \
#     "$PBMCS"  2>&1 |tee > featurecounts_rg.log 


# https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/bam


cd "$WD"/data


#bamtools split -in $PBMCS -tag CB
# crashes due to https://github.com/pezmaster31/bamtools/issues/135
samtools view $PBMCS | python split_script.py

