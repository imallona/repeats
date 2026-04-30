#!/usr/bin/env python
"""
Phase 2 prep: resolve SRA accessions per library for the repogle_2022 use case.

Reads selected_libraries.tsv (output of Phase 1) and queries SRA metadata
via pysradb to find the SRR run accessions for each KD6 library's mRNA
fastqs. The figshare library id (e.g. KD6_1_essential) maps to one or more
SRA records whose library_name follows the pattern

  KD6_seq{seq_run}_essential_mRNA_lane_{gem_group}_S{x}_L{lane}

The "lane_{gem_group}" token in the filename naming is the figshare gemgroup
id, not a sequencing lane. The actual sequencing lane is encoded in L00X.

Outputs:
  --out-yaml     YAML snippet usable as a second --configfile to the workflow
                 sc pipeline; populates dataset.samples.{library}.srrs.
  --out-tsv      selected_libraries.tsv augmented with srrs and n_srrs columns.
"""

import argparse
import csv
import re
import sys

import yaml

try:
    from pysradb import SRAweb
except ImportError as e:
    sys.stderr.write(f"pysradb not available: {e}\n")
    sys.exit(1)


def load_libraries_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fetch_sra_metadata(bioproject):
    db = SRAweb()
    return db.sra_metadata(bioproject, detailed=False)


def match_library_runs(df, gem_group, kd_label, source_pattern):
    pattern = re.compile(
        rf"^{kd_label}_seq\d+_essential_{source_pattern}_lane_{gem_group}_"
    )
    mask = df["library_name"].apply(
        lambda x: bool(pattern.match(x)) if isinstance(x, str) else False
    )
    return sorted(df.loc[mask, "run_accession"].dropna().unique().tolist())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bioproject", default="PRJNA831566")
    ap.add_argument("--libraries", required=True,
                    help="selected_libraries.tsv from Phase 1")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--workflow-config", required=True,
                    help="path to the workflow-side sc config yaml "
                         "(workflow/configs/repogle_kd6_sc.yaml); the merged "
                         "config emitted at --out-merged-yaml has the "
                         "dataset.samples block injected here")
    ap.add_argument("--out-merged-yaml", required=True,
                    help="path for the complete workflow config yaml with "
                         "samples merged in; this is the single --configfile "
                         "consumed by the workflow Phase 2 alignment")
    ap.add_argument("--kd-label", default="KD6",
                    help="library prefix used in SRA library_name (KD6 for "
                         "K562 essential, KD8 for K562 gwps)")
    ap.add_argument("--source-pattern", default="mRNA",
                    choices=["mRNA", "sgRNA"],
                    help="which read modality to resolve")
    args = ap.parse_args()

    libs = load_libraries_tsv(args.libraries)
    print(f"loaded {len(libs)} libraries from {args.libraries}", file=sys.stderr)

    print(f"querying SRA for bioproject {args.bioproject}", file=sys.stderr)
    df = fetch_sra_metadata(args.bioproject)
    print(f"  fetched {len(df)} SRA records", file=sys.stderr)

    samples = {}
    augmented = []
    for lib in libs:
        gg = lib["gem_group"]
        library_id = lib["library"]
        srrs = match_library_runs(df, gg, args.kd_label, args.source_pattern)
        samples[library_id] = {"srrs": srrs}
        out_lib = dict(lib)
        out_lib["srrs"] = ";".join(srrs)
        out_lib["n_srrs"] = len(srrs)
        augmented.append(out_lib)
        print(f"  {library_id} (gemgroup {gg}): {len(srrs)} SRRs", file=sys.stderr)

    if augmented:
        keys = list(augmented[0].keys())
        with open(args.out_tsv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t")
            w.writeheader()
            w.writerows(augmented)

    with open(args.workflow_config) as fh:
        merged = yaml.safe_load(fh)
    merged.setdefault("dataset", {})["samples"] = samples
    with open(args.out_merged_yaml, "w") as fh:
        yaml.safe_dump(merged, fh, sort_keys=False, default_flow_style=False)

    n_total = sum(len(s["srrs"]) for s in samples.values())
    print(f"wrote {args.out_tsv}", file=sys.stderr)
    print(f"wrote {args.out_merged_yaml} ({len(samples)} samples, {n_total} SRRs total)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
