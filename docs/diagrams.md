# workflow diagrams

Render to SVG with the mermaid CLI:

```
npm install -g @mermaid-js/mermaid-cli
mmdc -i diagrams.md -o diagrams.svg          # renders first diagram
mmdc -i diagrams.md -o diagrams/ --outputFormat svg  # one SVG per diagram
```

Or paste any block into https://mermaid.live for interactive editing and PNG/SVG export.

---

## overall pipeline

```mermaid
flowchart TD
    subgraph refs["reference preparation (reference.snmk)"]
        R1[genome FASTA\nGRCh38]
        R2[gene annotation\nEnsembl GTF]
        R3[repeat annotation\nRepeatMasker RMSK]
        R1 & R3 --> FA[feature FASTAs\nrepeats.fa\ngenic_repeats.fa\nintergenic_repeats.fa]
        R3 --> LM[locus maps\ntranscript_id / gene_id\nfamily_id / class_id]
        R1 & R2 & R3 --> IDX[aligner indices\nSTAR / Kallisto\nSalmon / Bowtie2]
        R2 & R3 --> FS[feature-set GTFs\ngenic_repeats.gtf\nintergenic_repeats.gtf]
    end

    subgraph sim["simulation (simulations.snmk)"]
        S1[draw expressed loci\nper cell from RMSK]
        S1 --> S2[sample reads from\ngenomic sequence]
        S2 --> FQ[FASTQ files\nSmartSeq2: one per cell\nChromium: R1 CB+UMI / R2 cDNA]
        S2 --> GT[ground_truth.tsv\ncell x locus x count]
    end

    subgraph aln["alignment and quantification"]
        A1[STARsolo\nstarsolo.snmk]
        A2[Kallisto + bustools\nkallisto.snmk]
        A3[Salmon alevin\nalevin.snmk]
        A4[Bowtie2 pseudo-genome\nbowtie2.snmk]
    end

    subgraph norm["normalization (normalize.snmk)"]
        N1[normalize_starsolo.py]
        N2[normalize_kallisto_*.py]
        N3[normalize_alevin_*.py]
        N4[count_pseudo_genome*.py]
    end

    subgraph eval["evaluation (evaluation.snmk)"]
        E1[evaluate.py\nglobal / per-cell / per-class]
        E2[aggregate_global_metrics]
        E3[evaluation_report.Rmd\nHTML report]
    end

    refs --> aln
    sim --> aln
    aln --> norm
    norm --> eval
    GT --> eval
    E1 --> E2 --> E3
```

---

## simulation design

```mermaid
flowchart LR
    subgraph input["inputs"]
        G[genome FASTA]
        R[repeat GTF\nRepeatMasker]
    end

    subgraph draw["per-cell locus sampling"]
        D1[draw n_expressed_per_cell loci\nfrom RMSK intervals]
        D2[assign read counts\nPoisson or fixed]
    end

    subgraph extract["sequence extraction"]
        E1[stream FASTA\none chromosome at a time]
        E2[extract subsequence\nfor each sampled locus]
    end

    subgraph ss2["SmartSeq2 output"]
        SS1[cell_001.fastq.gz\ncell_002.fastq.gz\n...\none file per cell]
    end

    subgraph chrom["Chromium output"]
        CR1[R1.fastq.gz\n16 nt CB + 12 nt UMI]
        CR2[R2.fastq.gz\ncDNA sequence]
        CR3[barcodes drawn from\n10x v3 whitelist\nunique UMI per cell+locus]
    end

    subgraph truth["ground truth"]
        GT[ground_truth.tsv\ncell_id / repeat_id\nfamily_id / class_id / true_count]
    end

    input --> draw --> extract
    extract --> ss2
    extract --> chrom
    draw --> truth
```

---

## bowtie2 pseudo-genome approach

```mermaid
flowchart TD
    subgraph build["index construction (once, shared across runs)"]
        GTF[repeat feature-set GTF]
        FA[genome FASTA]
        GTF & FA --> BED[BED of repeat loci]
        BED & FA --> PG["pseudo-genome FASTA\none entry per locus\ntranscript_id::chrom:start-end(strand)"]
        PG --> BT2IDX[bowtie2 index]
        GTF --> LM[locus_map.tsv\ntranscript_id / gene_id / family_id / class_id]
    end

    subgraph ss2path["SmartSeq2 path (per cell)"]
        FQ1[cell_NNN.fastq.gz]
        FQ1 --> AL1[bowtie2 -a\nall alignments]
        AL1 --> BAM1[cell BAM]
        BAM1 --> FILT[filter NH==1 for unique mode]
        FILT --> IDX1[samtools idxstats]
        IDX1 --> CPG[count_pseudo_genome.py\nlocus map lookup + granularity aggregation]
        CPG --> TSV1[counts_granularity_mode.tsv\nfeature x cell]
    end

    subgraph crpath["Chromium path (all cells)"]
        R2[R2.fastq.gz\ncDNA]
        R1[R1.fastq.gz\nCB + UMI]
        R2 --> AL2[bowtie2 -a\nall alignments]
        AL2 --> BAM2[aligned.bam]
        R1 --> TAG[tag_bam_chromium.py\nattach CB and UB tags]
        BAM2 --> TAG
        TAG --> TBAM[tagged.bam]
        TBAM --> DEDUP[umi_tools dedup --per-cell\nUMI deduplication]
        DEDUP --> DBAM[dedup.bam\nCB-tagged deduplicated]
        DBAM --> CPGC[count_pseudo_genome_chromium.py\nsingle-pass CB tag scan\nsparse accumulator]
        CPGC --> TSV2[counts_granularity.tsv\nfeature x cell]
    end

    BT2IDX --> AL1
    BT2IDX --> AL2
    LM --> CPG
    LM --> CPGC
```

---

## quantification granularity

```mermaid
flowchart LR
    subgraph lm["locus_map.tsv (4 columns, no header)"]
        COL["transcript_id | gene_id | family_id | class_id\nAluSz6::chr10:1000-1300+  |  AluSz6  |  Alu  |  SINE"]
    end

    subgraph gran["aggregation levels"]
        G0[locus\ntranscript_id\ne.g. AluSz6::chr10:1000-1300+]
        G1[gene_id\ne.g. AluSz6]
        G2[family_id\ne.g. Alu]
        G3[class_id\ne.g. SINE]
    end

    subgraph agg["aggregation operation"]
        A1[identity\nno summing]
        A2[sum counts across all loci\nwith same gene_id]
        A3[sum counts across all loci\nwith same family_id]
        A4[sum counts across all loci\nwith same class_id]
    end

    lm --> G0 --> A1
    lm --> G1 --> A2
    lm --> G2 --> A3
    lm --> G3 --> A4

    A1 & A2 & A3 & A4 --> OUT[feature x cell TSV\none per granularity]
```

---

## evaluation design

```mermaid
flowchart TD
    subgraph inputs["inputs"]
        OBS[observed counts\ncounts/aligner_fset_gran_mm.tsv\nfeature x cell]
        GT[ground_truth.tsv\ncell_id / repeat_id / true_count]
        LM2[locus_map.tsv\ndefines full feature universe\nfor specificity calculation]
        BM[snakemake benchmark .txt\nwall time / RSS / I/O]
    end

    subgraph agg2["ground truth aggregation"]
        AGG[aggregate to same granularity\nas observed counts]
    end

    subgraph metrics["evaluate.py metrics"]
        GM[global metrics\nPearson r, Spearman r, log1p RMSE\nprecision, recall, F1, Jaccard\nspecificity, compute resources]
        PC[per-cell metrics\nPearson r, Spearman r per cell]
        PF[per-class metrics\nall global metrics stratified\nby RepeatMasker class_id]
    end

    subgraph report["report"]
        SUM[summary_global_metrics.tsv\nall aligners concatenated]
        HTML[evaluation_report.html\nggplot2 bar charts\nprecision-recall scatter\nper-cell violin\nresource bar charts]
    end

    GT --> AGG --> GM & PC & PF
    OBS --> GM & PC & PF
    LM2 --> GM
    BM -.->|optional| GM
    GM --> SUM --> HTML
    PC --> HTML
    PF --> HTML
```
