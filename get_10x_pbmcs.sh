#!/bin/bas
##
## part of the repeats scrnaseq profiler fast prototype
##
## Izaskun Mallona
## 11th July 2019

WD=/home/Shared_s3it/imallona/repeats_sc
cd $WD

mkdir -p data
cd data
wget http://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/3.0.2/5k_pbmc_v3/5k_pbmc_v3_possorted_genome_bam.bam


# get hg19
