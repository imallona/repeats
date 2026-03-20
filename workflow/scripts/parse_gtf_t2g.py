#!/usr/bin/env python3
"""
Parse a GTF file and write a TSV mapping transcript_id to one or more
attribute values.

Usage examples
--------------
# 2-column t2g (genes GTF):
  parse_gtf_t2g.py --gtf genes.gtf --output genes_t2g.tsv --values gene_id

# 4-column locus_map (repeats GTF):
  parse_gtf_t2g.py --gtf repeats.gtf --output locus_map.tsv \
      --values gene_id family_id class_id

Works with any GTF that has transcript_id in the attribute field (column 9,
1-based).  Handles both Ensembl gene GTFs (gene_id / transcript_id) and
RepeatMasker-derived GTFs (transcript_id / gene_id / family_id / class_id).
"""
import argparse
import sys


def parse_attrs(attr_str):
    """Return dict of {key: value} from a GTF attribute string."""
    attrs = {}
    for part in attr_str.split(";"):
        part = part.strip()
        if not part:
            continue
        kv = part.split(" ", 1)
        if len(kv) == 2:
            attrs[kv[0]] = kv[1].strip('"')
    return attrs


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--feature", default="transcript",
                    help="GTF feature type to extract (column 3). Default: transcript")
    ap.add_argument("--key", default="transcript_id",
                    help="Attribute to use as the first (key) column. Default: transcript_id")
    ap.add_argument("--values", nargs="+", required=True,
                    help="Attribute(s) to use as value column(s). "
                         "E.g. --values gene_id   OR   --values gene_id family_id class_id")
    args = ap.parse_args()

    seen = set()
    with open(args.gtf) as fh, open(args.output, "w") as out:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if parts[2] != args.feature:
                continue
            attrs = parse_attrs(parts[8])
            key_val = attrs.get(args.key, "")
            if not key_val:
                continue
            row = (key_val,) + tuple(attrs.get(v, "") for v in args.values)
            if row not in seen:
                seen.add(row)
                out.write("\t".join(row) + "\n")


if __name__ == "__main__":
    main()
