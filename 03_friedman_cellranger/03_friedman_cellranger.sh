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
NTHREADS=20
MEM_GB=100
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

# we expected 5k cells each, as depicted
# https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6268/E-MTAB-6268.idf.txt
# but one of the samples gets around 18k cells, running it again specifying the expected
# number of cells

# ==> ./Day0Rep1_cellranger/outs/metrics_summary.csv <==
# Estimated Number of Cells,Mean Reads per Cell,Median Genes per Cell,Number of Reads,Valid Barcodes,Sequencing Saturation,Q30 Bases in Barcode,Q30 Bases in RNA Read,Q30 Bases in Sample Index,Q30 Bases in UMI,Reads Mapped to Genome,Reads Mapped Confidently to Genome,Reads Mapped Confidently to Intergenic Regions,Reads Mapped Confidently to Intronic Regions,Reads Mapped Confidently to Exonic Regions,Reads Mapped Confidently to Transcriptome,Reads Mapped Antisense to Gene,Fraction Reads in Cells,Total Genes Detected,Median UMI Counts per Cell
# "18,702","18,990","1,347","355,158,219",92.7%,46.9%,81.3%,81.2%,89.8%,84.2%,96.1%,93.6%,3.5%,9.9%,80.2%,76.1%,0.7%,80.3%,"22,382","3,382"

item=Day0Rep1
"$CELLRANGER" count --id="$item"_cellranger_5000_cells_expected \
              --fastqs="$DATA" \
              --transcriptome="$GRCh38" \
              --jobmode=local \
              --localcores=$NTHREADS \
              --sample="$item" \
              --localmem="$MEM_GB" \
              --expect-cells=5000 | tee "$item"_cellranger/run.log
