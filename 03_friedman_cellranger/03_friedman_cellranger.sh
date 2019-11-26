#!/bin/bash
##
## Runs a stardad cellranger count with Friedman's cardiac lineage differentiation
## Mouse genome
##
## run in portmac
## 
## Izaskun Mallona
## 25th Nov 2019


TAG="friedman_single_cell"
HOME=/home/imallona
DATA="/home/imallona/repeats_sc/data/""$TAG"
NTHREADS=30
MEM_GB=200
CELLRANGER="$HOME"/"soft/cellranger/cellranger-3.1.0/cellranger"
WD="$HOME"/"repeats_sc/data/friedman_single_cell"
GRCh38=$HOME"/indices/cellranger/refdata-cellranger-GRCh38-3.0.0"

mkdir -p $WD
cd $_


# https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/fastq-input

for item in Day0Rep1 Day0Rep2 \
                     Day15Rep1 Day15Rep2 \
                     Day2Rep1 Day2Rep2 \
                     Day30Rep1 Day30Rep2 \
                     Day5Rep1 Day5Rep2
do
    echo $item
    ln -s "$item"_I1.fastq.gz "$item"_S1_L001_I1_001.fastq.gz
    ln -s "$item"_R1.fastq.gz "$item"_S1_L001_R1_001.fastq.gz
    ln -s "$item"_R2.fastq.gz "$item"_S1_L001_R2_001.fastq.gz
    ln -s "$item"_R3.fastq.gz "$item"_S1_L001_R3_001.fastq.gz

    # mkdir -p "$item"_cellranger
    
    "$CELLRANGER" count --id="$item"_cellranger \
                  --fastqs="$DATA" \
                  --transcriptome="$GRCh38" \
                  --jobmode=local \
                  --localcores=$NTHREADS \
                  --sample="$item" \
                   --localmem="$MEM_GB" | tee "$item"_cellranger/run.log

done 
