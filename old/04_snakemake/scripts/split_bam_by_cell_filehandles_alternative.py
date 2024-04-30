#! /usr/bin/env python3
#
# Modified from https://raw.githubusercontent.com/crazyhottommy/scATACutils/master/python/split_scATAC_bam_by_cell.py
#
# Gets a BAM file with the reads from the first NUMBEROFCELLS_MAX different CB (cell barcodes)
# found when iterating the bam records

import pysam
import csv
import argparse
import os.path
import sys

NUMCONCURRENT = 1000

parser = argparse.ArgumentParser()
parser.add_argument("bam", help="Required. the FULL path to the 10x scATAC bam file generated \
    by cellranger-atac count and sorted by CB")
parser.add_argument("-prefix", help="Optional, the prefix of the output bam, default is barcode.bam")
parser.add_argument("-outdir", help="Optional, the output directory for the splitted bams, default is current dir")
args = parser.parse_args()


if os.path.exists(args.bam):
    pass
else:
    print("bam not found")
    sys.exit(1)

if os.path.isdir(args.outdir):
    pass
else:
    try:
        os.mkdir(args.outdir)
    except OSError:
        print("can not create directory {}".format(args.outdir))

fin = pysam.AlignmentFile(args.bam, "rb")

## dict with keys are barcode, values are outbam
fouts_dict = {}

for read in fin:
    tags = read.tags
    CB_list = [ x for x in tags if x[0] == "CB"]
    if CB_list:
        cell_barcode = CB_list[0][1]
        if args.prefix:
            fout_name = args.prefix + "_" + cell_barcode + ".bam"
        else:
            fout_name = cell_barcode + ".bam"
        # print(len(fouts_dict))
        # if the number of filehandles is too high, close one
        if len(fouts_dict) >= NUMCONCURRENT and cell_barcode not in fouts_dict:
            fout_dict[fout_dict.keys()[1]].close()
            
        if len(fouts_dict) < NUMCONCURRENT:
            if cell_barcode not in fouts_dict:
                if args.outdir:
                    fouts_dict[cell_barcode] = pysam.AlignmentFile(os.path.join(args.outdir,fout_name), "wb", template = fin)
                else:
                    fouts_dict[cell_barcode] = pysam.AlignmentFile(fout_name, "wb", template = fin)
            if cell_barcode in fouts_dict:
                fouts_dict[cell_barcode].write(read)
    else: 
        continue
    

## do not forget to close the files
fin.close()
for fout in fouts_dict.values():
    fout.close()
