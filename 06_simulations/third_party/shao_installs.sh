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

# to be run

# ~/soft/TE_expression_in_scRNAseq/scripts/zUMIs_TE_modified/zUMIs-master-TE.sh -y zUMIs_sample.yaml -d path_to_modified_zUMIs_package

# ~/soft/TE_expression_in_scRNAseq/scripts/scRNA_UMI_counting.R -f Sample_1_featurecount.sam.txt -r mm10_pc_exon_te_tx.saf -n Sample_1 -s 50 -m 1 -p 10
