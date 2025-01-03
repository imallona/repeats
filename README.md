# Status

- currently (2024): active development, see [TODO](TODO.md)
- [old version v0.1](https://github.com/imallona/repeats/releases/tag/v0.1)

# Summary

We aimed to quantify conventional features (GENCODE) as well as repetitive elements (RepeatMasker) from multiple single-cell RNA-seq datasets using several analytical approaches, and starting from raw data, and looked for association with known cell metadata (e.g. cell types).

# Running

```
snakemake --use-conda --cores 10 -p
```

# v0.1 release

## Data

We reanalyzed the datasets listed below. We are very grateful to the original data producers and encourage checking their data release licenses and citing their publication(s).

Numbers point to the specific snakemake workflow.

Data are available at [Zenodo](https://doi.org/10.5281/zenodo.4584956) as generated by `04_snakemake/export_run_for_upload.R`.

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

# Contact

izaskun dot mallona at gmail dot com or izaskun dot mallona at mls dot uzh dot ch

# License

GNU General Public License (GPL)

# Dates

- Started Thu Jul 11 14:21:11 CEST 2019
- Funded 15 Nov 2019
- Re-enabled Mon Nov 18 14:18:09 CET 2019
- Re-enabled Mon Oct 14 01:40:21 PM CEST 2024

