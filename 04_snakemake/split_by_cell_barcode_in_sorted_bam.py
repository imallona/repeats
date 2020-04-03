#! /usr/bin/env python3
#
import pysam
import csv
import argparse
import os.path
import sys

parser = argparse.ArgumentParser()
parser.add_argument("bam", help="BAM file with sorted RG tags")
parser.add_argument("-outdir", help="Tthe output directory for the splitted bams, default is current dir")

args = parser.parse_args()

if os.path.exists(args.bam):
    pass
else:
    print("bam not found")
    sys.exit(1)

fin = pysam.AlignmentFile(args.bam, "rb")

if os.path.isdir(args.outdir):
    pass
else:
    try:
        os.mkdir(args.outdir)
    except OSError:
        print("failed to make directory {}".format(args.outdir))

prev_cell_barcode = 'empty'

for read in fin:
    tags = read.tags
    RG_list = [ x for x in tags if x[0] == "RG"]
    if RG_list:
        cell_barcode = RG_list[0][1]

        if prev_cell_barcode is 'empty':
            fout_name = cell_barcode + ".bam"
            fout_name = os.path.join(args.outdir, fout_name)
            
            fout = pysam.AlignmentFile(fout_name, "wb", template = fin)            
        elif cell_barcode == prev_cell_barcode:
            fout.write(read)
        elif cell_barcode != prev_cell_barcode:
            prev_fout.close()
            fout_name = cell_barcode + ".bam"
            fout_name = os.path.join(args.outdir, fout_name)
            fout = pysam.AlignmentFile(fout_name, "wb", template = fin)            
            fout.write(read)
        
        prev_fout = fout
        prev_cell_barcode = cell_barcode
    else: 
        continue
    
fin.close()
fout.close()
