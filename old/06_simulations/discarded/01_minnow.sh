#!/bin/bash
#
# 26 Apr 2021
# https://github.com/COMBINE-lab/minnow#TLDR

export MINNOW="~/soft/minnow/build/src/minnow"
KMER_LENGTH=101

# build salmon index (skipped, we have them)

# mkdir data

# wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.pc_transcripts.fa.gz
# gunzip data/gencode.vM25.pc_transcripts.fa.gz

# wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz
# gunzip data/gencode.vM25.annotation.gtf.gz

# wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz
# gunzip data/GRCm38.primary_assembly.genome.fa.gz

# grep ">" GRCm38.primary_assembly.genome.fa | cut -d ">" -f 2 | cut -d " " -f 1 > GRCm38.primary_assembly.genome.chrnames.txt

# cd ..

# salmon index \
# -t gencode.vM25.pc_transcripts.fa \
# -i gencode.vM25.annotation.sidx --gencode -p 128 \
# -d GRCm38.primary_assembly.genome.chrnames.txt



## Step 0.2 -- Run Alevin on a mouse dataset (skipped, we have alevin runs)

# cd data
# wget https://sra-pub-src-1.s3.amazonaws.com/SRR6459157/AdultMouse_Rep3_possorted_genome_bam.bam.1

# mv AdultMouse_Rep3_possorted_genome_bam.bam.1 AdultMouse_Rep3_possorted_genome_bam.bam

# bamtofastq --reads-per-fastq=500000000 AdultMouse_Rep3_possorted_genome_bam.bam FASTQtmp

# mv FASTQtmp/Ad-Ms-Total-Sorted_20k_count_MissingLibrary_1_HK2GNBBXX/bamtofastq_S1_L006_I1_001.fastq.gz AdultMouseRep3_S1_L001_I1_001.fastq.gz

# mv FASTQtmp/Ad-Ms-Total-Sorted_20k_count_MissingLibrary_1_HK2GNBBXX/bamtofastq_S1_L006_R1_001.fastq.gz AdultMouseRep3_S1_L001_R1_001.fastq.gz

# mv FASTQtmp/Ad-Ms-Total-Sorted_20k_count_MissingLibrary_1_HK2GNBBXX/bamtofastq_S1_L006_R2_001.fastq.gz AdultMouseRep3_S1_L001_R2_001.fastq.gz
# cd ..


# awk -F "\t" '$3 == "transcript" { print $9 }' data/gencode.vM25.annotation.gtf | tr -d ";\"" | awk '{print $4"\t"$2}' > data/gencode.vM25.annotation.tx2gene.tsv

# salmon alevin -l ISR -i gencode.vM25.annotation.sidx \
# -1 data/AdultMouseRep3_S1_L001_R1_001.fastq.gz \
# -2 data/AdultMouseRep3_S1_L001_R2_001.fastq.gz \
# -o alevin_out -p 36 \
# --tgMap data/gencode.vM25.annotation.tx2gene.tsv \
# --chromium \
# --dumpFeatures --expectCells 1850 \
# --dumpBfh


# but the dumpBfh we don't have, run it with a modified version, rule `map_salmon_repeats_chromium`


# (snakemake) [imallona@imlsportmacquarie 04_snakemake]$ snakemake -s 02_repeats_pbmc5k_chromium.snmk  -n -p ~/repeats_sc/runs/5k_pbmc_v3/alevin/repeats/alevin/quants_mat.gz 


# mv ~/repeats_sc/runs/5k_pbmc_v3/alevin/repeats/alevin/quants_mat.gz{,.old}
# snakemake -s 02_repeats_pbmc5k_chromium.snmk  -n -p ~/repeats_sc/runs/5k_pbmc_v3/alevin/repeats/alevin/quants_mat.gz

# mv ~/repeats_sc/runs/5k_pbmc_v3/alevin/repeats/alevin/quants_mat.gz.old ~/repeats_sc/runs/5k_pbmc_v3/alevin/repeats/alevin/quants_mat.gz


mkdir -p /home/imallona/repeats_sc/runs/5k_pbmc_v3/minnow

cd /home/imallona/repeats_sc/runs/5k_pbmc_v3/minnow

(~/soft/salmon/salmon-1.1.0_linux_x86_64/bin/salmon alevin     -l ISR     -1 /home/imallona/repeats_sc/data/5k_pbmc_v3/5k_pbmc_v3_fastqs/5k_pbmc_v3_S1_L001_R1_001.fastq.gz /home/imallona/repeats_sc/data/5k_pbmc_v3/5k_pbmc_v3_fastqs/5k_pbmc_v3_S1_L002_R1_001.fastq.gz /home/imallona/repeats_sc/data/5k_pbmc_v3/5k_pbmc_v3_fastqs/5k_pbmc_v3_S1_L003_R1_001.fastq.gz /home/imallona/repeats_sc/data/5k_pbmc_v3/5k_pbmc_v3_fastqs/5k_pbmc_v3_S1_L004_R1_001.fastq.gz     -2 /home/imallona/repeats_sc/data/5k_pbmc_v3/5k_pbmc_v3_fastqs/5k_pbmc_v3_S1_L001_R2_001.fastq.gz /home/imallona/repeats_sc/data/5k_pbmc_v3/5k_pbmc_v3_fastqs/5k_pbmc_v3_S1_L002_R2_001.fastq.gz /home/imallona/repeats_sc/data/5k_pbmc_v3/5k_pbmc_v3_fastqs/5k_pbmc_v3_S1_L003_R2_001.fastq.gz /home/imallona/repeats_sc/data/5k_pbmc_v3/5k_pbmc_v3_fastqs/5k_pbmc_v3_S1_L004_R2_001.fastq.gz  \
                                                    --chromiumV3     \
                                                    -i /home/imallona/repeats_sc/indices/salmon/GRCh38/repeats_salmon   \
                                                    -p 40 \
                                                    -o /home/imallona/repeats_sc/runs/5k_pbmc_v3/minnow/repeats    \
                                                    --dumpBfh \
                                                    --tgMap /home/imallona/repeats_sc/indices/salmon/GRCh38/repeats_salmon/txp2gene.tsv ) 2> /home/imallona/repeats_sc/runs/5k_pbmc_v3/minnow/run_salmon_repeats_chromium.log
.gz



## 1.1 De-Bruijn graph construction using minnow index

# getting an index from the repeatome itself!
# -k kmer length, important to keep it fixed for the simulation
# -f, --filt-size for twopaco index (?)
~/soft/minnow/build/src/minnow index -r /home/imallona/repeats_sc/annotation/repeatome_from_mm10_rmsk_TE.gtf.fa.gz \
           -k "$KMER_LENGTH" \
           -f 20 --tmpdir ~/tmp \
           -p 10 \
           -o /home/imallona/repeats_sc/indices/minnow/repeats


# till here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## minnow estimate

# crashes
~/soft/minnow/build/src/minnow estimate \
                               -o ~/repeats_sc/runs/5k_pbmc_v3/minnow \
                               -r ~/repeats_sc/annotation/repeatome_from_mm10_rmsk_TE.gtf.fa.gz \
                               --g2t ~/repeats_sc/indices/salmon/mm10/repeats_salmon/txp2gene.tsv \
                               --bfh /home/imallona/repeats_sc/runs/5k_pbmc_v3/minnow/repeats/alevin/bfh.txt



gunzip ~/repeats_sc/annotation/repeatome_from_mm10_rmsk_TE.gtf.fa.gz

# crashes too
~/soft/minnow/build/src/minnow estimate \
                               -o ~/repeats_sc/runs/5k_pbmc_v3/minnow \
                               -r ~/repeats_sc/annotation/repeatome_from_mm10_rmsk_TE.gtf.fa \
                               --g2t ~/repeats_sc/indices/salmon/mm10/repeats_salmon/txp2gene.tsv \
                               --bfh /home/imallona/repeats_sc/runs/5k_pbmc_v3/minnow/repeats/alevin/bfh.txt


pigz -p 10 ~/repeats_sc/annotation/repeatome_from_mm10_rmsk_TE.gtf.fa


crashes at this level, maybe try to simulate it myself?


