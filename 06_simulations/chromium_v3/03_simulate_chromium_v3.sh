#!/bin/bash


POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -l|--length)
    LENGTH="$2"
    shift # past argument
    shift # past value
    ;;
    -n|--nthreads)
    NTHREADS="$2"
    shift # past argument
    shift # past value
    ;;
    -q|--qual_symbol)
    QUAL_SYMBOL="$2"
    shift # past argument
    shift # past value
    ;;
    -g|--gtf)
    GTF="$2"
    shift    
    shift 
    ;;
    --genome_fasta)
    GENOME_FASTA="$2"
    shift    
    shift 
    ;;
    --bedtools_bin)
    BEDTOOLS_BIN="$2"
    shift    
    shift 
    ;;
    --i1_origin)
    I1_ORIGIN="$2"
    shift    
    shift 
    ;;
    --r1_origin)
    R1_ORIGIN="$2"
    shift    
    shift 
    ;;
    --r2_sim_origin)
    R2_SIM_ORIGIN="$2"
    shift    
    shift 
    ;;
    --i1_sim)
    I1_SIM="$2"
    shift    
    shift 
    ;;
    --r1_sim)
    R1_SIM="$2"
    shift    
    shift 
    ;;
    --truth_count_table)
    TRUTH_COUNT_TABLE="$2"
    shift    
    shift 
    ;;
    --cb)
    CB="$2"
    shift    
    shift 
    ;;
    --umi)
    UMI="$2"
    shift    
    shift 
    ;;
    --num_cells)
    NUM_CELLS="$2"
    shift    
    shift 
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

if [[ -n $1 ]]; then
    echo "Last line of file specified as non-opt/last argument:"
    tail -1 "$1"
fi




# LENGTH=91 # beware, this will generate 92 nt long fastqs
# NTHREADS=30
# QUAL_SYMBOL="F"
# # get rmsk coordinates
# GTF="$HOME/repeats_sc/annotation/GRCh38_rmsk_TE.gtf.gz"

# GENOME_FASTA="$HOME/repeats_sc/annotation/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

# BEDTOOLS="$HOME/soft/bedtools/bedtools-2.29.2/bin/bedtools"


# I1_ORIGIN="$HOME/repeats_sc/data/5k_pbmc_v3/5k_pbmc_v3_fastqs/5k_pbmc_v3_S1_L002_I1_001.fastq.gz"
# R1_ORIGIN="$HOME/repeats_sc/data/5k_pbmc_v3/5k_pbmc_v3_fastqs/5k_pbmc_v3_S1_L002_R1_001.fastq.gz"
# R2_SIM_ORIGIN="$HOME/repeats_sc/data/sim_5k_pbmc_v3/sim_5k_pbmc_v3_S1_L002_R2_001.fastq"


# I1_SIM="$HOME/repeats_sc/data/sim_5k_pbmc_v3/sim_5k_pbmc_v3_S1_L002_I1_001.fastq"
# R1_SIM="$HOME/repeats_sc/data/sim_5k_pbmc_v3/sim_5k_pbmc_v3_S1_L002_R1_001.fastq"
# TRUTH_COUNT_TABLE="$HOME/repeats_sc/data/sim_5k_pbmc_v3/truth.tsv"

# CB=16
# UMI=12

# NUM_CELLS=1000

source 02_simulate_ground_truth_fastqs.sh

mkdir -p $HOME/repeats_sc/data/sim_5k_pbmc_v3

## function from the script above
touch "$R2_SIM_ORIGIN"
echo simulate_cdna $LENGTH $QUAL_SYMBOL $GTF $GENOME_FASTA $NTHREADS $BEDTOOLS_BIN $R2_SIM_ORIGIN


simulate_cdna $LENGTH $QUAL_SYMBOL $GTF $GENOME_FASTA $NTHREADS $BEDTOOLS_BIN $R2_SIM_ORIGIN

gunzip $HOME/repeats_sc/data/sim_5k_pbmc_v3/*fastq.gz

# postprocessing
/usr/local/R/R-4.0.5/bin/Rscript simulate_chromium.R \
                                 $R1_ORIGIN \
                                 $R2_SIM_ORIGIN \
                                 $I1_SIM \
                                 $R1_SIM \
                                 $NTHREADS \
                                 $CB \
                                 $UMI \
                                 $NUM_CELLS \
                                 $TRUTH_COUNT_TABLE

pigz -p $NTHREADS $HOME/repeats_sc/data/sim_5k_pbmc_v3/*fastq
pigz -p $NTHREADS $HOME/repeats_sc/data/sim_5k_pbmc_v3/truth.tsv


# # add 5 mutations to each fastq record - test
## caution untested!
# /usr/local/R/R-4.0.5/bin/Rscript add_mutations_to_fastq.R \
#                                  $R2_SIMORIGIN \
#                                  5 \
#                                  test_5.fastq
