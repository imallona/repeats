#!/bin/bash

# target=04_repeats_pbmc_10k_v3_chromium.snmk
target=09_pbmc8k_chromium.snmk
cp "$target"{,.backup}

sed -i 's/BASE/config["base"]/g' "$target"
sed -i 's/RUN_NAME/config["run_name"]/g' "$target"
sed -i 's/CHEMISTRY/config["chemistry"]/g' "$target"
sed -i 's/LOCAL_MEM_GB/config["params"]["local_mem_gb"]/g' "$target"
sed -i 's/NTHREADS/config["params"]["nthreads"]/g' "$target"
sed -i 's/GENOME_SHORT/config["genome_short"]/g' "$target"
sed -i 's/GENOME_URL/config["genome_url"]/g' "$target"
sed -i 's/GENOME/config["genome"]/g' "$target"
sed -i 's/RUN_NAME/config["run_name"]/g' "$target"
sed -i 's/FEATURECOUNTS_RSCRIPT/config["dependencies"]["featurecounts_parsing"]/g' "$target"
sed -i 's/REP_GTF_URL/config["rep_gtf_url"]/g' "$target"
sed -i 's/GENES_GTF_URL/config["genes_gtf_url"]/g' "$target"


# sed -i 's/BEDOPS/config["software"]["bedops"]/g' "$target"
