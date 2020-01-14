# Summary

Understanding the determinants of cell identity is a large challenge in biology, especially with the uptake of new technologies. Cell identity is  intimately linked to transcriptional regulation. Overall protein coding gene expression is well accepted as a proxy to cell identity, whilst the repetitive elements expression contribution to cell phenotypes is largely unknown.

Taking advantage of the single-cell RNA-Seq (scRNA-seq) technology we propose a methodological development to profile repeat expression and to cluster cells according to repeats transcriptional landscapes. Largely abundant in eukaryotic genomes, repeats are expressed and show regulatory signatures; studies carried out in bulk RNAseq indicate that repeat transcription is linked to cell identity.

To further dissect the repeat expression landscape in single cells, we will compare cell-to-cell variability arising from repeat expression to that of protein-coding genes in adult brain cells (with strongly defined cell types), in a model for cardiac differentiation (transcriptional trajectories) and in cancer (dysregulated repeat expression).

Methodologically, we will address the sequencing, mappability, quantification and normalization challenges and will implement an open, robust and flexible computational biology workflow.

As a proof of concept to test the project feasibility we have run a prototype on a droplet-based (10x) dataset of peripheral mononuclear blood cells (PBMCs), known by the shallow sequencing and 3' RNA enrichment. We detected a noticeable amount of repeat expression and variability across cell identities.

From the impact point of view, we will open a window into the overlooked transcriptomics of repetitive DNA in single-cells, and will release a tool to dissect any present or future dataset in any species.

# Files

by exec order

`00_prototype`
- get_10x_pbmcs.sh    gets 10x pbmc dataset
- get_gtf.sh          gets grh38 rmsk table
- feature_counts.sh   counts rmsk instances from 10x pbmc table

`01_descriptive`
- `01_descriptive.Rmd` gets the Spark application preliminary data

`02_data_retrieval` retrieves the datasets from the Spark application (mouse)

# Remotes and branch strategy

No-fast-forward `develop` to `master` (with `hotfix`) branching strategy

- origin	git@bitbucket.org:imallona/repeats_sc.git (fetch)
- origin	git@bitbucket.org:imallona/repeats_sc.git (push)
- github	https://github.com/imallona/repeats_sc.git (fetch)
- github	https://github.com/imallona/repeats_sc.git (push)

# For benchmark

- TEToolkit http://hammelllab.labsites.cshl.edu/software/
- RepeatScout http://bix.ucsd.edu/repeatscout/

# Dates

- Started Thu Jul 11 14:21:11 CEST 2019
- Funded 15 Nov 2019
- Re-enabled Mon Nov 18 14:18:09 CET 2019
