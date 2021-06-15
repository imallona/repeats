#!/usr/bin/env R

library(data.table)

gtf_fn <- '~/repeats_sc/annotation/GRCh38_rmsk_TE.gtf.gz'

gtf <- fread(gtf_fn, colClasses = c(rep('NULL', 8), 'character'))

reg <- 'gene_id "(.+)"; transcript_id "(.+)"; family_id "(.+)"; class_id "(.+)";'

gtf$id <- gsub(pattern = reg, replacement = '\\2', x = gtf$V9)
gtf$name <- gsub(pattern = reg, replacement = '\\1', x = gtf$V9)
gtf$family <- gsub(pattern = reg, replacement = '\\3', x = gtf$V9)
gtf$class <- gsub(pattern = reg, replacement = '\\4', x = gtf$V9)

gtf <- as.data.frame(gtf[,c('id', 'name', 'family', 'class')])
## rownames(gtf) <- gtf$id
## gtf <- gtf[,-1]

head(gtf)
# read the files and check
fn <- '/home/imallona/repeats_sc/runs/sim_5k_pbmc_v3/count_repeats_loci_on_cellranger_standard/multimappers/cellranger_standard:0:1:cb:TTTGACTCAGATACTC-1.bam.featureCounts'

foo <- as.data.frame(fread(fn, sep = '\t'))
table(foo$V2)

## > table(foo$V2)
##              Assigned  Unassigned_Ambiguity Unassigned_NoFeatures 
##                  3822                   948                    18 
##   Unassigned_Unmapped 
## 24

fn <- '/home/imallona/repeats_sc/runs/sim_5k_pbmc_v3/count_repeats_loci_on_cellranger_standard/unique_reads/cellranger_standard:0:1:cb:TTTGACTCAGATACTC-1.bam.featureCounts'

foo <- as.data.frame(fread(fn, sep = '\t'))
table(foo$V2)
## > table(foo$V2)

##                Assigned    Unassigned_Ambiguity Unassigned_MultiMapping 
##                    3488                     891                     393 
##   Unassigned_NoFeatures     Unassigned_Unmapped 
##                      16                      24 
## >


## fn <- '~/repeats_sc/runs/sim_5k_pbmc_v3/test_count_repeats_loci_on_cellranger_standard/multimappers/cellranger_standard:0:1:cb:TTTGATCCAACCACAT-1.bam.featureCounts'
fn <- '~/repeats_sc/runs/sim_5k_pbmc_v3/test_count_repeats_loci_on_cellranger_standard/multimappers/cellranger_standard:0:1:cb:TTTCGATCAATTTCGG-1.bam.featureCounts'

## samtools view "/home/imallona/repeats_sc/runs/sim_5k_pbmc_v3/count_repeats_on_cellranger_standard/split/cellranger_standard:0:1:cb:TTTCGATCAATTTCGG-1.bam" | less

foo <- as.data.frame(fread(fn, sep = '\t'))
head(foo)

## 1 ACRO1;ACRO1_dup63;acro;Satellite;::GL000214.1:137003-137094 Assigned  1
## 2              Alu;Alu_dup1412;Alu;SINE;::6:47128686-47128777 Assigned  2
## 3              Alu;Alu_dup1849;Alu;SINE;::8:14452472-14452563 Assigned  1
## 4            Alu;Alu_dup2000;Alu;SINE;::8:140609862-140609953 Assigned  2
## 5           Alu;Alu_dup2993;Alu;SINE;::12:103286944-103287035 Assigned  2
## 6               Alu;Alu_dup462;Alu;SINE;::2:72550908-72550999 Assigned  1
##                           V4
## 1                ACRO1_dup63
## 2 Alu_dup1412,MLT1A0_dup8044
## 3                Alu_dup1849
## 4 Alu_dup2000,AluJo_dup35071
## 5  Alu_dup2993,L1PA4_dup9747
## 6                 Alu_dup462



## > table(foo$V2)

##              Assigned Unassigned_NoFeatures   Unassigned_Unmapped 
## 4698                    17                    34

table(foo$V3)

## foo$mapped_id <- gsub(pattern = '(.+);(.+);(.+);(.+);.*',
##                    replacement = '\\2', x = foo$V1)

## foo$mapped_name <- gsub(pattern = '(.+);(.+);(.+);(.+);.*',
##                       replacement = '\\1', x = foo$V1)

## foo$mapped_family <- gsub(pattern = '(.+);(.+);(.+);(.+);.*',
##                       replacement = '\\3', x = foo$V1)
## foo$mapped_class <- gsub(pattern = '(.+);(.+);(.+);(.+);.*',
##                       replacement = '\\4', x = foo$V1)


foo <- merge(foo, gtf, by.x= 'V4', by.y = 'id', all.x = TRUE)
colnames(foo) <- c('mapped.id', 'read_name', 'status', 'count', 'mapped.name', 'mapped.family',
                   'mapped.class')

foo$sim.id <- gsub(pattern = '(.+);(.+);(.+);(.+);.*', replacement = '\\2', x = foo$read_name)
foo <- merge(foo, gtf, by.x= 'sim.id', by.y = 'id', all.x = TRUE)
head(foo)

colnames(foo)[9:11] <- c('sim.name', 'sim.family', 'sim.class') 

table(foo$mapped.family == foo$sim.family, useNA = 'always')
table(foo$mapped.name == foo$sim.name, useNA = 'always')
table(foo$mapped.class == foo$sim.class, useNA = 'always')
table(foo$status)



## what???
head(foo)
gtf_full <- as.data.frame(fread(gtf_fn))
gtf_full[grepl('Alu_dup2578', gtf_full$V9),]
gtf_full[grepl('L1ME3A_dup9744', gtf_full$V9),]

## 1       7SK_dup431                 7SK_dup431
## 2 ALR/Alpha_dup987           ALR/Alpha_dup987
## 3       Alu_dup108                 Alu_dup108
## 4      Alu_dup2578 Alu_dup2578,L1ME3A_dup9744
## 5      Alu_dup3068 Alu_dup3068,AluJo_dup53239
## 6      Alu_dup3102                Alu_dup3102

## > gtf_full[grepl('Alu_dup2578', gtf_full$V9),]
##        V1        V2   V3       V4       V5  V6 V7 V8
## 496229 10 hg38_rmsk exon 77131324 77131374 311  +  .
##                                                                                   V9
## 496229 gene_id "Alu"; transcript_id "Alu_dup2578"; family_id "Alu"; class_id "SINE";
## > gtf_full[grepl('L1ME3A_dup9744', gtf_full$V9),]
##        V1        V2   V3       V4       V5  V6 V7 V8
## 496230 10 hg38_rmsk exon 77131411 77131682 918  +  .
##                                                                                        V9
## 496230 gene_id "L1ME3A"; transcript_id "L1ME3A_dup9744"; family_id "L1"; class_id "LINE";




### let's spawn items

