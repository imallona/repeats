#!/bin/bash

find ~/repeats_sc/runs -name "*FALSE*FALSE.rds" | grep -v 'cobra\|aris\|markers'

# /home/imallona/repeats_sc/runs/ding_seqwell/ding_seqwell_pmbc_seqwell_regress_nCount_FALSE_nFeature_FALSE.rds
# /home/imallona/repeats_sc/runs/ding_dropseq/ding_dropseq_pmbc_dropseq_regress_nCount_FALSE_nFeature_FALSE.rds
# /home/imallona/repeats_sc/runs/SRR10974767/SRR10974767_chromium_regress_nCount_FALSE_nFeature_FALSE.rds
# /home/imallona/repeats_sc/runs/pbmc8k/pbmc8k_pmbc_chromium_regress_nCount_FALSE_nFeature_FALSE.rds
# /home/imallona/repeats_sc/runs/pbmc8k/pbmc8k_pmbc_alevin_regress_nCount_FALSE_nFeature_FALSE.rds
# /home/imallona/repeats_sc/runs/frozen_pbmc_donor_a/summarize_cellranger_run_regress_nCount_FALSE_nFeature_FALSE.rds
# /home/imallona/repeats_sc/runs/frozen_pbmc_donor_a/frozen_pbmc_donor_a_pmbc_cellranger_regress_nCount_FALSE_nFeature_FALSE.rds
# /home/imallona/repeats_sc/runs/frozen_pbmc_donor_a/frozen_pbmc_donor_a_pmbc_chromium_regress_nCount_FALSE_nFeature_FALSE.rds
# /home/imallona/repeats_sc/runs/5k_pbmc_v3/summarize_cellranger_run_regress_nCount_FALSE_nFeature_FALSE.rds
# /home/imallona/repeats_sc/runs/5k_pbmc_v3/5k_pbmc_v3_pmbc_chromium_regress_nCount_FALSE_nFeature_FALSE.rds
# /home/imallona/repeats_sc/runs/5k_pbmc_v3/5k_pbmc_v3_pmbc_alevin_regress_nCount_FALSE_nFeature_FALSE.rds
# /home/imallona/repeats_sc/runs/SRR8847571/SRR8847571_chromium_regress_nCount_FALSE_nFeature_FALSE.rds
# /home/imallona/repeats_sc/runs/SRR10974769/summarize_cellranger_run_regress_nCount_FALSE_nFeature_FALSE.rds
# /home/imallona/repeats_sc/runs/SRR10974769/SRR10974769_pmbc_cellranger_regress_nCount_FALSE_nFeature_FALSE.rds
# /home/imallona/repeats_sc/runs/SRR10974769/SRR10974769_chromium_regress_nCount_FALSE_nFeature_FALSE.rds
# /home/imallona/repeats_sc/runs/pbmcs_smartseq2_ding/pbmcs_smartseq2_ding_pmbc_smartseq2_regress_nCount_FALSE_nFeature_FALSE.rds
# /home/imallona/repeats_sc/runs/zheng_truth/zheng_truth_pbmc_zheng_truth_regress_nCount_FALSE_nFeature_FALSE.rds
# /home/imallona/repeats_sc/runs/ding_celseq2/ding_celseq2_pmbc_celseq2_regress_nCount_FALSE_nFeature_FALSE.rds
# /home/imallona/repeats_sc/runs/pbmc_10k_v3/pbmc_10k_v3_pmbc_alevin_regress_nCount_FALSE_nFeature_FALSE.rds
# /home/imallona/repeats_sc/runs/pbmc_10k_v3/pbmc_10k_v3_pmbc_chromium_regress_nCount_FALSE_nFeature_FALSE.rds

data=/home/imallona/repeats_sc/runs/
exp=/home/imallona/src/repeats_sc/04_snakemake/export_run_for_upload.R
mkdir -p ~/zenodo
cd $_

Rscript "$exp" --input  "$data"/ding_seqwell/ding_seqwell_pmbc_seqwell_regress_nCount_FALSE_nFeature_FALSE.rds --output ding_seqwell.tar
Rscript "$exp" --input    "$data"/ding_dropseq/ding_dropseq_pmbc_dropseq_regress_nCount_FALSE_nFeature_FALSE.rds --output ding_dropseq.tar
Rscript "$exp" --input               "$data"/SRR10974767/SRR10974767_chromium_regress_nCount_FALSE_nFeature_FALSE.rds --output SRR10974767.tar
Rscript "$exp" --input               "$data"/pbmc8k/pbmc8k_pmbc_chromium_regress_nCount_FALSE_nFeature_FALSE.rds --output pbmc8k.tar
Rscript "$exp" --input               "$data"/frozen_pbmc_donor_a/frozen_pbmc_donor_a_pmbc_chromium_regress_nCount_FALSE_nFeature_FALSE.rds \
        --output frozen_pbmc_donor_a.tar
Rscript "$exp" --input               "$data"/5k_pbmc_v3/5k_pbmc_v3_pmbc_chromium_regress_nCount_FALSE_nFeature_FALSE.rds --output 5k_pbmc_v3.tar
Rscript "$exp" --input               "$data"/SRR8847571/SRR8847571_chromium_regress_nCount_FALSE_nFeature_FALSE.rds --output SRR8847571.tar
Rscript "$exp" --input               "$data"/SRR10974769/SRR10974769_chromium_regress_nCount_FALSE_nFeature_FALSE.rds --output SRR10974769.tar
Rscript "$exp" --input               "$data"/pbmcs_smartseq2_ding/pbmcs_smartseq2_ding_pmbc_smartseq2_regress_nCount_FALSE_nFeature_FALSE.rds --output pbmcs_smartseq2_ding.tar
Rscript "$exp" --input               "$data"/zheng_truth/zheng_truth_pbmc_zheng_truth_regress_nCount_FALSE_nFeature_FALSE.rds --output zheng_truth.tar
Rscript "$exp" --input               "$data"/ding_celseq2/ding_celseq2_pmbc_celseq2_regress_nCount_FALSE_nFeature_FALSE.rds --output ding_celseq2.tar
Rscript "$exp" --input               "$data"/pbmc_10k_v3/pbmc_10k_v3_pmbc_chromium_regress_nCount_FALSE_nFeature_FALSE.rds --output pbmc_10k_v3.tar
