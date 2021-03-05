#!/usr/bin/env R
##
## Exports Seurat object as compressedplain text files for zenodo
##   (sparse or not)
##
## 5 March 2021
## Izaskun Mallona
## GPL

# README file

readme <- "
# Purpose

Quantification outputs and embeddings from [repeats_sc](https://bitbucket.org/imallona/repeats_sc/src/master/). Files report conventional features (GENCODE) as well as repetitive elements (rmsk) from several sources and using multiple analytical approaches.

# Files

Each `dataset` and `flavour` (see definitions below) has the following files:

1. `md5sum`: the original Seurat object md5sum
2. `counts.mtx.gz`: raw counts (features: rows, cells: columns), matrixmarket
3. `norm_counts.mtx.gz`: normalized counts (features: rows, cells: columns), matrixmarket
4. `cell_meta.csv.gz`: cell metadata, comma-separated
5. `feature_meta.csv.gz`: feature metadata, comma-separated
6. `pca_cell_embeddings.csv.gz`, PCA embedding, comma-separated
7. `pca_feature.loadings.csv.gz`, PCA loadings, comma-separated
8. `tsne_cell_embeddings.csv.gz`, t-SNE embedding, comma-separated
9. `umap_cell_embeddings.csv.gz`, UMAP embedding, comma-separated

# Datasets

We reanalyzed the datasets listed below. We are very grateful to the original data producers and encourage checking their data release licenses and citing their publication(s).

Numbers point to the specific snakemake workflow.

- `5k_pbmc_v3`, `02`, [10x](https://cf.10xgenomics.com/samples/cell-exp/3.0.2/5k_pbmc_v3/5k_pbmc_v3_web_summary.html)
- `brain`, `17`, [SRA](https://sra-pub-src-1.s3.amazonaws.com/SRR6854141/10X50_3.bam.1)
- `castration`, `18`, [SRA](https://www.ncbi.nlm.nih.gov/sra?term=SRP256199)
- `ding_celseq2`, `08`, [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132044)
- `ding_dropseq`, `14`, [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132044)
- `ding_seqwell`, `15`, [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132044)
- `frozen_pbmc_donor_a`, `05`, [10x](http://cf.10xgenomics.com/samples/cell-exp/1.1.0/frozen_pbmc_donor_a/frozen_pbmc_donor_a_web_summary.html)
- `GSE121861_kumar_mouse`, `13`, [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121861)
- `pbmc_10k_v3`, `04`, [10x](https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_web_summary.html)
- `pbmc8k`, `09`, [10x](https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc8k)
- `pbmcs_smartseq2_ding`, `06`, [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132044)
- `SRR10974767`
- `SRR10974769`, `01`, [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4286760)
- `SRR8847571`, `10`, [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4286760)
- `SRR9307700`, `11`, [SRA](https://www.ncbi.nlm.nih.gov/sra?term=SRX7639834)
- `zheng_truth`, `16`, [10x](https://support.10xgenomics.com/single-cell-gene-expression/datasets)

# Flavours

Whether the analyte are genes or repeats or both; whether the mapper was alevin, star, or bowtie; weather the multimappers are counted or not; and whether the repeats might overlap genes or not.

Please notice some datasets are not suitable for some analytical strategies; for instance, we didn't run `alevin` on Smart-seq2 data.

- `both_alevin`
- `both_multi`
- `both_unique`
- `control_repeats_only_multi`
- `control_repeats_only_unique`
- `genes_alevin`
- `genes_cellranger`
- `repeats_alevin`
- `repeats_bowtie_multi`
- `repeats_bowtie_unique`
- `repeats_cellranger`
- `repeats_star_multi`
- `repeats_star_unique`

# Repository

- [Bitbucket](https://bitbucket.org/imallona/repeats_sc/src/master/)

# Preprint

- In preparation

# Contact

izaskun.mallona at gmail.com

# Funding

- [Swiss National Science Foundation](http://p3.snf.ch/project-190824)
Institut für Molekulare Biologie (IMLS) Universität Zürich
"

suppressPackageStartupMessages({
    library("argparse")
    library(Seurat)
    library(Matrix)
    library(R.utils)
})

parser <- ArgumentParser()

parser$add_argument("-i", "--input",type  =  'character',
    help = "List of Seurat objects (RDS)")

parser$add_argument("-o", "--output", type = 'character', help = "Output tarball (.tar)")

args <- parser$parse_args()

## print(args)
if (is.null(args$input) | is.null(args$output))
    stop('Bad arguments')
## args <- list(input = '~/repeats_sc/runs/pbmc8k/pbmc8k_pmbc_chromium_regress_nCount_TRUE_nFeature_TRUE.rds',
##              output = '~/tmp/pbmc8k.tgz')
## print(args)

d <- readRDS(args$input)

outdir <- file.path(tempdir(), gsub('.tar', '', basename(args$output)),
                    gsub('.tar', '', basename(args$output)))

## print (outdir)

## input seurat object md5sum
dir.create(outdir, recursive = TRUE)
system(sprintf('md5sum %s > %s.md5sum', args$input, file.path(outdir)))


for (flavour in names(d)) {

    ## used to have line breaks - removing them here
    for (item in grep('orig', colnames(d[[flavour]]@meta.data)))
        d[[flavour]]@meta.data[,item] <- gsub('\n', '', d[[flavour]]@meta.data[,item])

        
    path <- file.path(outdir, flavour)
    dir.create(path, recursive = TRUE)
    
    ## counts
    writeMM(obj = Matrix(GetAssayData(d[[flavour]], 'counts'), sparse = TRUE),
            file = file.path(path, 'counts.mtx.gz'))
    
    gzip(file.path(path, 'counts.mtx.gz'), overwrite = TRUE)

    ## normalized counts
    writeMM(obj = Matrix(GetAssayData(d[[flavour]], 'data'), sparse = TRUE),
            file = file.path(path, 'norm_counts.mtx.gz'))
    
    gzip(file.path(path, 'norm_counts.mtx.gz'), overwrite = TRUE)

    ## cell meta
    gzmeta <- gzfile(file.path(path, "cell_meta.csv.gz"), "w")
    write.csv(d[[flavour]]@meta.data, gzmeta, row.names = FALSE)
    close(gzmeta)

    ## feature metadata (just the featurename)
    gzmeta <- gzfile(file.path(path, "feature_meta.csv.gz"), "w")
    write.csv(data.frame(feature = row.names(d[[flavour]]@meta.data)), gzmeta,
              row.names = FALSE)
    close(gzmeta)

    ## pca
    ## embeddings
    gz <- gzfile(file.path(path, "pca_cell_embeddings.csv.gz"), "w")
    write.csv(d[[flavour]]@reductions$pca@cell.embeddings, gz, row.names = FALSE)
    close(gz)

    ## feature.loadings
    gz <- gzfile(file.path(path, "pca_feature.loadings.csv.gz"), "w")
    write.csv(d[[flavour]]@reductions$pca@feature.loadings, gz, row.names = FALSE)
    close(gz)
    
    ## tsne
    gz <- gzfile(file.path(path, "tsne_cell_embeddings.csv.gz"), "w")
    write.csv(d[[flavour]]@reductions$tsne@cell.embeddings, gz, row.names = FALSE)
    close(gz)
    
    ## umap
    gz <- gzfile(file.path(path, "umap_cell_embeddings.csv.gz"), "w")
    write.csv(d[[flavour]]@reductions$umap@cell.embeddings, gz, row.names = FALSE)
    close(gz)

}

## tar(args$output, files = outdir,
##     compression = c("none"),
##     tar = Sys.getenv("tar"),
##     extra_flags = sprintf("-C %s", outdir))

system(sprintf('tar cvf %s -C %s ..', args$output, outdir))
