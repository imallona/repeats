#!/bin/bash
#
# https://github.com/wanqingshao/TE_expression_in_scRNAseq
#
# 3 May 2021


# install stringtie2

cd ~/soft

git clone https://github.com/skovaka/stringtie2

cd stringtie2
make release

# taco

mkdir -p ~/soft/taco

cd $_

wget https://github.com/tacorna/taco/releases/download/v0.7.3/taco-v0.7.3.Linux_x86_64.tar.gz

tar xzvf taco-v0.7.3.Linux_x86_64.tar.gz


cd ~/soft

git clone https://github.com/wanqingshao/TE_expression_in_scRNAseq



# Rpackages install

BiocManager::install(c('data.table',
                       'GenomicRanges',
                       'stringidst',
'ggplot2',
'ggpubr',
'ggrepel',
'lattice',
'methods',
'reshape2',
'rtracklayer',
'yaml',
'data.table',
'GenomicRanges',
'ggplot2',
'ggpubr',
'ggrepel',
'magrittr',
'optparse',
'parallel'))


###


# First, generate SAF files

# GENES=~/repeats_sc/annotation/gencode.v33.primary_assembly.annotation.gtf.gz
# # REPS=~/repeats_sc/annotation/GRCh38_rmsk_TE.gtf.gz
# Rscript=/usr/local/R/R-4.1.0/bin/Rscript
# # remove the chr starts for $GENES, since $REPS doesn't have them

# zcat $GENES | sed 's/^chr//g'  > ~/tmp/gencode_v33_nochr.gtf

wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/rmsk.txt.gz -O ~/tmp/rmsk.txt.gz
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/refGene.txt.gz -O ~/tmp/refGene.txt.gz


mysql  --user=genome --host=genome-mysql.cse.ucsc.edu -A -D hg38 -e 'select size from chromInfo' -B -N > ~/tmp/chrlengths.txt 


# zcat ~/tmp/rmsk.txt.gz > ~/tmp/rmsk.txt
# zcat ~/tmp/refGene.txt.gz > ~/tmp/refGene.txt

$Rscript ~/soft/TE_expression_in_scRNAseq/scripts/export_saf_files.r \
         -g ~/tmp/refGene.txt.gz \
         -r ~/tmp/rmsk.txt.gz \
         -l ~/tmp/chrlengths.txt \
         -n shao_GRCh38 \
         -m chrM -p 20

mkdir -p ~/repeats_sc/annotation/shao
mv *saf ~/repeats_sc/annotation/shao



cat << EOF > /home/imallona/src/repeats_sc/06_simulations/third_party/zumis2.yaml
###########################################
#Welcome to zUMIs
#below, please fill the mandatory inputs
#We expect full paths for all files.
###########################################

#define a project name that will be used to name output files
project: sim_5k_zumis

#Sequencing File Inputs:
#For each input file, make one list object & define path and barcode ranges
#base definition vocabulary: BC(n) UMI(n) cDNA(n).
#Barcode range definition needs to account for all ranges. You can give several comma-separated ranges for BC & UMI sequences, eg. BC(1-6,20-26)
#you can specify between 1 and 4 input files
sequence_files:
  file1:
    name: /home/imallona/repeats_sc/data/sim_5k_pbmc_v3/sim_5k_pbmc_v3_S1_L002_R2_001.fastq.gz  #path to first file
    base_definition:
      - BC(1-16) #example: BC(1-6)
      - UMI(17-28) #example: UMI(7-16)
  file2:
    name: /home/imallona/repeats_sc/data/sim_5k_pbmc_v3/sim_5k_pbmc_v3_S1_L002_R1_001.fastq.gz  #path to second file
    base_definition:
      - cDNA(1-92) #example: cDNA(1-50)

#reference genome setup
reference:
  STAR_index: /home/imallona/repeats_sc/indices/star/GRCh38/transcriptome/
  SAF_file: /home/imallona/repeats_sc/annotation/shao/shao_GRCh38_te.saf
  # GTF_file: /home/imallona/repeats_sc/annotation/gencode.v33.primary_assembly.annotation.gtf.gz
  additional_files: #Optional parameter. It is possible to give additional reference sequences here, eg ERCC.fa
  additional_STAR_params: #Optional parameter. you may add custom mapping parameters to STAR here
  allow_multimapping: yes #yes or no for multimapping option
#output directory
out_dir: /home/imallona/repeats_sc/runs/sim_5k_pbmc_v3/zumis # directory

###########################################
#below, you may optionally change default parameters
###########################################

#number of processors to use
num_threads: 20
mem_limit: 100 #Memory limit in Gigabytes, null meaning unlimited RAM usage.

#barcode & UMI filtering options
#number of bases under the base quality cutoff that should be filtered out.
#Phred score base-cutoff for quality control.
filter_cutoffs:
  BC_filter:
    num_bases: 1
    phred: 20
  UMI_filter:
    num_bases: 1
    phred: 20

#Options for Barcode handling
#You can give either number of top barcodes to use or give an annotation of cell barcodes.
#If you leave both barcode_num and barcode_file empty, zUMIs will perform automatic cell barcode selection for you!
barcodes:
  barcode_num: null
  barcode_file: null
  automatic: yes #Give yes/no to this option. If the cell barcodes should be detected automatically. If the barcode file is given in combination with automatic barcode detection, the list of given barcodes will be used as whitelist.
  BarcodeBinning: 1 #Hamming distance binning of close cell barcode sequences.
  nReadsperCell: 100 #Keep only the cell barcodes with atleast n number of reads.

#Options related to counting of reads towards expression profiles
counting_opts:
  strand: 1 #Is the library stranded? 0 = unstranded, 1 = positively stranded, 2 = negatively stranded
  countMultiMappingReads: yes #Do you want to count the multiple aligned reads
  fracOverlap: 0.1 #percentage reads to overlap with features
  twoPass: yes #perform basic STAR twoPass mapping


#Start zUMIs from stage. Possible TEXT(Filtering, Mapping, Counting, Summarising). Default: Filtering.
which_Stage: Filtering

#define dependencies program paths
samtools_exec: samtools #samtools executable
Rscript_exec: /usr/local/R/R-4.1.0/bin/Rscript # Rscript #Rscript executable
STAR_exec: /home/imallona/soft/star/STAR-2.7.3a/bin/Linux_x86_64/STAR #STAR executable
pigz_exec: pigz #pigz executable

#below, fqfilter will add a read_layout flag defining SE or PE
zUMIs_directory: /home/imallona/soft/TE_expression_in_scRNAseq/scripts/zUMIs_TE_modified/
read_layout: SE
EOF


mkdir -p ~/repeats_sc/runs/sim_5k_pbmc_v3/zumis

~/soft/TE_expression_in_scRNAseq/scripts/zUMIs_TE_modified/zUMIs-master-TE.sh \
    -y /home/imallona/src/repeats_sc/06_simulations/third_party/zumis2.yaml \
    -d /home/imallona/soft/TE_expression_in_scRNAseq/scripts/zUMIs_TE_modified

"$Rscript" ~/soft/TE_expression_in_scRNAseq/scripts/scRNA_UMI_counting.R \
    -f ~/repeats_sc/runs/sim_5k_pbmc_v3/zumis/sim_5k_zumis_featurecount.sam.txt \
    -r ~/repeats_sc/annotation/shao/shao_GRCh38_te.saf \
    -n sim_5k_zumis \
    -s 50 \
    -m 1 \
    -p 30

# The final output for is Sample_1_count_dt.txt. This file contains five columns: "feature" is the feature name (genes and TE transcripts), "uniq" quantifies the number of reads that are mapped to single features, "EM_distri" quantifies the number of reads at each feature after distributing reads that are mapped to multiple features using EM algorithm, "even_distri" quantifies the number of reads at each feature after evenly distributing reads that are mapped to multiple features, "RG" specifies the cell barcode.

# [imallona@imlsportmacquarie third_party]$ head sim_5k_zumis_count_dt.txt
# feature uniq    EM_distri       even_distri     RG
# te_1000999      0       0.0989610178227491      0.333333333333333       AAAAAAAAAAAAAAAA
# te_1001307      1       1       1       AAAAAAAAAAAAAAAA
# te_1004288      1       1       1       AAAAAAAAAAAAAAAA
# te_1028134      0       0.672645739910314       0.5     AAAAAAAAAAAAAAAA
# te_1028195      0       0.327354260089686       0.5     AAAAAAAAAAAAAAAA
# te_1029295      0       0.185290658697604       0.25    AAAAAAAAAAAAAAAA
# te_1085630      0       0.123223395205381       0.166666666666667       AAAAAAAAAAAAAAAA
# te_1103335      1       1       1       AAAAAAAAAAAAAAAA
# te_1115112      2       2       2       AAAAAAAAAAAAAAAA

# and ~/repeats_sc/annotation/shao/shao_GRCh38_te.saf has the start-end (?)

mv sim_5k_zumis_count_dt.txt ~/repeats_sc/runs/sim_5k_pbmc_v3/zumis/
