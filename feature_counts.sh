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
samtools view $PBMCS | /usr/bin/python3 bam_cb_split.py 

# mkdir -p bamtools
# cd $_
# ln -s ../data/5*bam

# ulimit -n 200000
# bamtools split -in $PBMCS -tag CB


# trying differently
# https://github.com/crazyhottommy/scATACutils/blob/master/python/split_scATAC_bam_by_cell.py

# requires a virtenv with pysam

cd ~/virtenvs
virtualenv -p python3 split_bam
# cd $_
 # cd split_bam
source bin/activate
pip -V
pip install pysam
# pip install csv

cd $WD/data

ulimit -n 50000
python3 split_bam_by_cell.py -prefix foo -outdir . $PBMCS

deactivate

# then samtobam and featurecounts
# # @todo implement
