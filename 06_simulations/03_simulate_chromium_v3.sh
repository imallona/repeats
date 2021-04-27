#!/bin/bash



LENGTH=92
NTHREADS=20
QUAL_SYMBOL="F"
# get rmsk coordinates
GTF="$HOME/repeats_sc/annotation/mm10_rmsk_TE.gtf.gz"

GENOME_FASTA="$HOME/repeats_sc/annotation/Mus_musculus.GRCm38.dna.primary_assembly.fa"

BEDTOOLS="$HOME/soft/bedtools/bedtools-2.29.2/bin/bedtools"


I1_ORIGIN="$HOME/repeats_sc/data/5k_pbmc_v3/5k_pbmc_v3_fastqs/5k_pbmc_v3_S1_L002_I1_001.fastq.gz"
R1_ORIGIN="$HOME/repeats_sc/data/5k_pbmc_v3/5k_pbmc_v3_fastqs/5k_pbmc_v3_S1_L002_R1_001.fastq.gz"
R2_SIM_ORIGIN="$HOME/repeats_sc/data/sim/5k_pbmc_v3_sim_S1_L002_R2_001.fastq"


I1_SIM="$HOME/repeats_sc/data/sim/5k_pbmc_v3_sim_S1_L002_I1_001.fastq"
R1_SIM="$HOME/repeats_sc/data/sim/5k_pbmc_v3_sim_S1_L002_R1_001.fastq"

CB=16
UMI=12

NUM_CELLS=1000

source 02_simulate_ground_truth_fastqs.sh

## function from the script above
touch "$R2_SIM_ORIGIN"
simulate_cdna $LENGTH $QUAL_SYMBOL $GTF $GENOME_FASTA $NTHREADS $BEDTOOLS $R2_SIM_ORIGIN

# postprocessing
/usr/local/R/R-4.0.5/bin/Rscript simulate_chromium.R


