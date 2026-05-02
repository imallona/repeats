#!/usr/bin/env python3
"""
Parse one or more STAR Log.final.out files and emit a tidy TSV summary.

Each input log path is expected to be at
  <starsolo_base>/<sample_id>/<mode>_<feature_set>/Log.final.out
so that sample_id, mode, and feature_set can be derived from the path.

Output columns:
  sample_id, multimapper_mode, feature_set, log_path,
  n_input_reads, avg_input_length, avg_mapped_length,
  uniquely_mapped_pct, multi_mapped_pct, too_many_loci_pct,
  unmapped_too_short_pct, unmapped_other_pct,
  mismatch_rate_pct, deletion_rate_pct, insertion_rate_pct,
  deletion_avg_length, insertion_avg_length

Percentages are kept in their reported scale (e.g. "5.32" not 0.0532).
Missing fields become empty cells. Designed to be called once per run with a
list of logs and aggregated by snakemake/R.
"""

import argparse
import os
import re
import sys


KEY_MAP = {
    "Number of input reads": "n_input_reads",
    "Average input read length": "avg_input_length",
    "Average mapped length": "avg_mapped_length",
    "Uniquely mapped reads %": "uniquely_mapped_pct",
    "% of reads mapped to multiple loci": "multi_mapped_pct",
    "% of reads mapped to too many loci": "too_many_loci_pct",
    "% of reads unmapped: too short": "unmapped_too_short_pct",
    "% of reads unmapped: other": "unmapped_other_pct",
    "Mismatch rate per base, %": "mismatch_rate_pct",
    "Deletion rate per base": "deletion_rate_pct",
    "Insertion rate per base": "insertion_rate_pct",
    "Deletion average length": "deletion_avg_length",
    "Insertion average length": "insertion_avg_length",
}

PCT_SUFFIX = re.compile(r"\s*%\s*$")


def parse_one(path):
    out = {v: "" for v in KEY_MAP.values()}
    with open(path) as fh:
        for line in fh:
            if "|" not in line:
                continue
            key, value = [s.strip() for s in line.split("|", 1)]
            target = KEY_MAP.get(key)
            if target is None:
                continue
            value = PCT_SUFFIX.sub("", value)
            out[target] = value
    return out


def derive_path_fields(log_path):
    """Recover (sample_id, multimapper_mode, feature_set) from the log path.

    Expected layout:
      .../starsolo/<sample_id>/<mode>_<feature_set>/Log.final.out
    Falls back to empty strings for non-matching layouts so callers can still
    consume the row.
    """
    parent = os.path.basename(os.path.dirname(log_path))
    grandparent = os.path.basename(
        os.path.dirname(os.path.dirname(log_path)))
    mode = ""
    feature_set = ""
    if "_" in parent:
        head, _, tail = parent.partition("_")
        if head in ("unique", "multi"):
            mode = head
            feature_set = tail
    return grandparent, mode, feature_set


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--log", action="append", required=True,
                    help="One or more Log.final.out paths")
    ap.add_argument("--out", required=True, help="Output TSV")
    args = ap.parse_args()

    columns = [
        "sample_id", "multimapper_mode", "feature_set", "log_path",
        *KEY_MAP.values(),
    ]

    with open(args.out, "w") as out:
        out.write("\t".join(columns) + "\n")
        for log_path in args.log:
            sample_id, mode, feature_set = derive_path_fields(log_path)
            try:
                row = parse_one(log_path)
            except FileNotFoundError:
                print(f"missing log: {log_path}", file=sys.stderr)
                continue
            row.update({
                "sample_id": sample_id,
                "multimapper_mode": mode,
                "feature_set": feature_set,
                "log_path": log_path,
            })
            out.write("\t".join(row[c] for c in columns) + "\n")


if __name__ == "__main__":
    main()
