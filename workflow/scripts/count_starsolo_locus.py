#!/usr/bin/env python3
"""
Count reads per (cell barcode, repeat locus) from a STARsolo BAM.

STARsolo aligns reads to the full genome. This script overlaps each alignment
with repeat locus intervals (parsed from the feature FASTA headers) to assign it
to a locus_id (= transcript_id, e.g. AluSz6_dup1).

Locus assignment: a read is assigned to a locus if the alignment start position
falls within the locus interval. For reads overlapping multiple loci, only the
first (smallest genomic start) matching locus is counted.

Interval coordinates are loaded from the feature FASTA whose headers encode
transcript_id::chrom:start-end(strand). This is a permanent index artifact,
unlike the repeats.gtf which is a temp file. The locus_id stored in the output
is the transcript_id (bare, no coordinates), matching locus_map col 0 and the
count matrices produced by bowtie2/kallisto/alevin.

Modes:
  smartseq2 - Cell ID is the prefix of the read name before _r{n}_,
              set by simulate_reads.py. All reads per (cell, locus) are counted.
              No CB:Z tag is present in SmartSeq2 STARsolo BAMs.
  chromium  - CB:Z = cell barcode, UB:Z = UMI. UMI deduplication is performed
              per (CB, locus) by counting distinct UMIs per cell.

Scalability:
  Chromium: BAM is sorted by CB tag into a temporary file. Reads are streamed
  one cell at a time, so peak memory is O(loci_per_cell x umis_per_locus).
  SmartSeq2: BAM is read directly (sorted by coordinate). All cells are
  accumulated simultaneously; with ~20 cells this is negligible memory.

Output: feature x cell TSV (rows=locus_ids, cols=cell barcodes, first col=feature_id).
"""

import argparse
import bisect
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict


def parse_locus_id(locus_id):
    """
    Parse a locus_id string of the form gene_id::chrom:start-end(strand).
    Returns (chrom, start_0, end_0, strand) or None if malformed.
    start_0 is 0-based, end_0 is 0-based exclusive (= GTF end field value).
    """
    try:
        gene_part, coords = locus_id.split("::", 1)
        colon_idx = coords.rfind(":")
        chrom = coords[:colon_idx]
        rest = coords[colon_idx + 1:]
        dash_idx = rest.index("-")
        start_0 = int(rest[:dash_idx])
        paren_idx = rest.index("(")
        end_0 = int(rest[dash_idx + 1:paren_idx])
        strand = rest[paren_idx + 1:rest.index(")")]
        return chrom, start_0, end_0, strand
    except (ValueError, IndexError):
        return None


def load_intervals(fasta_path):
    """
    Parse a feature FASTA file whose headers encode coordinates as
    transcript_id::chrom:start-end(strand) and build per-chrom sorted
    interval lists for fast overlap queries.

    locus_id is the transcript_id (part before ::), which is the canonical
    per-locus identifier used in locus_map col 0 and count matrices.

    Returns:
      chrom_starts: {chrom: sorted list of start positions}
      chrom_intervals: {chrom: list of (start, end, locus_id) sorted by start}
    """
    chrom_intervals = defaultdict(list)
    with open(fasta_path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            header = line[1:].rstrip("\n").split()[0]
            if "::" not in header:
                continue
            locus_id = header.split("::", 1)[0]
            coords = parse_locus_id(header)
            if coords is None:
                continue
            chrom, start_0, end_0, strand = coords
            chrom_intervals[chrom].append((start_0, end_0, locus_id))

    chrom_starts = {}
    for chrom in chrom_intervals:
        chrom_intervals[chrom].sort(key=lambda x: x[0])
        seen = set()
        deduped = []
        for entry in chrom_intervals[chrom]:
            if entry not in seen:
                seen.add(entry)
                deduped.append(entry)
        chrom_intervals[chrom] = deduped
        chrom_starts[chrom] = [iv[0] for iv in chrom_intervals[chrom]]

    return chrom_starts, dict(chrom_intervals)


def find_locus(chrom, pos, chrom_starts, chrom_intervals):
    """
    Find the first locus whose interval contains pos (0-based).
    Returns locus_id string or None.
    """
    intervals = chrom_intervals.get(chrom)
    if not intervals:
        return None
    starts = chrom_starts[chrom]
    idx = bisect.bisect_right(starts, pos) - 1
    for i in range(max(0, idx), min(len(intervals), idx + 5)):
        s, e, lid = intervals[i]
        if s > pos:
            break
        if s <= pos < e:
            return lid
    return None


def sort_bam_by_cb(bam_path, n_threads):
    """Sort BAM by CB tag into a tempfile. Returns tempfile path (caller must unlink)."""
    tmpf = tempfile.NamedTemporaryFile(suffix=".bam", delete=False)
    sorted_path = tmpf.name
    tmpf.close()
    print(f"Sorting BAM by CB tag -> {sorted_path}", file=sys.stderr)
    result = subprocess.run(
        ["samtools", "sort", "-t", "CB", "-@", str(n_threads), "-o", sorted_path, bam_path],
        stderr=subprocess.PIPE
    )
    if result.returncode != 0:
        os.unlink(sorted_path)
        sys.stderr.write(result.stderr.decode())
        sys.exit(f"samtools sort failed (exit {result.returncode})")
    return sorted_path


def process_smartseq2(bam_path, chrom_starts, chrom_intervals, multimapper_mode):
    """
    Process SmartSeq2 BAM: cell_id is the prefix of QNAME before _r{n}_.
    Reads are sorted by coordinate (interleaved across cells) so we accumulate
    all cells simultaneously.
    """
    counts = defaultdict(lambda: defaultdict(int))
    all_cbs = set()
    all_loci = set()
    n_reads = 0
    n_assigned = 0

    proc = subprocess.Popen(
        ["samtools", "view", bam_path],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    try:
        for raw in proc.stdout:
            line = raw.decode()
            if line.startswith("@"):
                continue
            fields = line.split("\t")
            if len(fields) < 11:
                continue
            flag = int(fields[1])
            if flag & 4 or flag & 2048:
                continue
            n_reads += 1

            nh = None
            for tag in fields[11:]:
                if tag.startswith("NH:i:"):
                    nh = int(tag[5:].rstrip("\n"))
                    break
            if multimapper_mode == "unique" and nh is not None and nh > 1:
                continue

            # cell_id is everything before _r{digits}_ in the QNAME
            cb = re.sub(r"_r\d+_.*$", "", fields[0])

            chrom = fields[2]
            pos_0 = int(fields[3]) - 1
            locus_id = find_locus(chrom, pos_0, chrom_starts, chrom_intervals)
            if locus_id is None:
                continue

            counts[cb][locus_id] += 1
            all_cbs.add(cb)
            all_loci.add(locus_id)
            n_assigned += 1
    finally:
        proc.kill()
        proc.wait()

    print(f"  {n_reads} reads processed, {n_assigned} assigned to loci", file=sys.stderr)
    print(f"  {len(all_cbs)} cells, {len(all_loci)} loci", file=sys.stderr)
    return dict(counts), all_cbs, all_loci


def process_chromium(bam_path, chrom_starts, chrom_intervals, multimapper_mode, n_threads):
    """
    Process Chromium BAM: sort by CB tag, stream one cell at a time to cap memory.
    UMI deduplication per (CB, locus).
    """
    sorted_bam = sort_bam_by_cb(bam_path, n_threads)

    counts = defaultdict(lambda: defaultdict(int))
    all_cbs = set()
    all_loci = set()
    n_reads = 0
    n_assigned = 0

    current_cb = None
    cell_data = {}

    proc = subprocess.Popen(
        ["samtools", "view", sorted_bam],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    try:
        for raw in proc.stdout:
            line = raw.decode()
            if line.startswith("@"):
                continue
            fields = line.split("\t")
            if len(fields) < 11:
                continue
            flag = int(fields[1])
            if flag & 4 or flag & 2048:
                continue
            n_reads += 1

            cb = None
            ub = None
            nh = None
            for tag in fields[11:]:
                if tag.startswith("CB:Z:"):
                    cb = tag[5:].rstrip("\n")
                elif tag.startswith("UB:Z:"):
                    ub = tag[5:].rstrip("\n")
                elif tag.startswith("NH:i:"):
                    nh = int(tag[5:].rstrip("\n"))

            if cb is None:
                continue
            if multimapper_mode == "unique" and nh is not None and nh > 1:
                continue

            locus_id = find_locus(fields[2], int(fields[3]) - 1, chrom_starts, chrom_intervals)
            if locus_id is None:
                continue

            if cb != current_cb:
                if current_cb is not None:
                    for lid, umis in cell_data.items():
                        counts[current_cb][lid] = len(umis)
                    all_cbs.add(current_cb)
                current_cb = cb
                cell_data = defaultdict(set)

            cell_data[locus_id].add(ub if ub else str(n_reads))
            all_loci.add(locus_id)
            n_assigned += 1
    finally:
        proc.kill()
        proc.wait()

    if current_cb is not None:
        for lid, umis in cell_data.items():
            counts[current_cb][lid] = len(umis)
        all_cbs.add(current_cb)

    os.unlink(sorted_bam)

    print(f"  {n_reads} reads processed, {n_assigned} assigned to loci", file=sys.stderr)
    print(f"  {len(all_cbs)} cells, {len(all_loci)} loci", file=sys.stderr)
    return dict(counts), all_cbs, all_loci


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bam", required=True,
                    help="STARsolo Aligned.sortedByCoord.out.bam")
    ap.add_argument("--fasta", required=True,
                    help="Feature FASTA with headers transcript_id::chrom:start-end(strand)")
    ap.add_argument("--mode", required=True, choices=["smartseq2", "chromium"])
    ap.add_argument("--multimapper-mode", default="unique",
                    choices=["unique", "multi"],
                    help="unique: NH=1 reads only; multi: all aligned reads")
    ap.add_argument("--threads", type=int, default=1,
                    help="Threads for samtools sort (Chromium only; default: 1)")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    print(f"Loading intervals from {args.fasta}", file=sys.stderr)
    chrom_starts, chrom_intervals = load_intervals(args.fasta)
    n_loci = sum(len(v) for v in chrom_intervals.values())
    print(f"  {n_loci} intervals on {len(chrom_intervals)} chromosomes", file=sys.stderr)

    if args.mode == "smartseq2":
        counts, all_cbs, all_loci = process_smartseq2(
            args.bam, chrom_starts, chrom_intervals, args.multimapper_mode)
    else:
        counts, all_cbs, all_loci = process_chromium(
            args.bam, chrom_starts, chrom_intervals, args.multimapper_mode, args.threads)

    sorted_loci = sorted(all_loci)
    sorted_cbs = sorted(all_cbs)

    with open(args.output, "w") as fh:
        fh.write("feature_id\t" + "\t".join(sorted_cbs) + "\n")
        for locus_id in sorted_loci:
            row = [str(counts.get(cb, {}).get(locus_id, 0)) for cb in sorted_cbs]
            fh.write(locus_id + "\t" + "\t".join(row) + "\n")

    print(f"wrote {len(sorted_loci)} loci x {len(sorted_cbs)} cells to {args.output}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
