#!/usr/bin/env python3
"""Split a GTF file into N chunks by gene_id, for parallel featureCounts."""
import argparse, math, os

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--gtf', required=True)
    ap.add_argument('--n-chunks', type=int, required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--prefix', required=True, help='Output filename prefix, e.g. repeats_feat_chunk')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    lines_by_gene = {}
    gene_order = []
    with open(args.gtf) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 9:
                continue
            gid = ''
            for attr in parts[8].split(';'):
                attr = attr.strip()
                if attr.startswith('gene_id'):
                    gid = attr[len('gene_id'):].strip().strip('"')
                    break
            if gid not in lines_by_gene:
                lines_by_gene[gid] = []
                gene_order.append(gid)
            lines_by_gene[gid].append(line)

    n = len(gene_order)
    csize = max(1, math.ceil(n / args.n_chunks))
    gid_to_chunk = {gid: min(i // csize, args.n_chunks - 1)
                    for i, gid in enumerate(gene_order)}

    out_paths = [os.path.join(args.outdir, f'{args.prefix}_{i}.gtf')
                 for i in range(args.n_chunks)]
    handles = [open(p, 'w') for p in out_paths]
    for gid, glines in lines_by_gene.items():
        idx = gid_to_chunk[gid]
        for gl in glines:
            handles[idx].write(gl)
    for h in handles:
        h.close()

if __name__ == '__main__':
    main()
