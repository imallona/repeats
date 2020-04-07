#! /usr/bin/env python3
#
import pysam
import csv
import argparse
import os.path
import sys
import gzip

parser = argparse.ArgumentParser()
parser.add_argument("bam", help="BAM file with sorted RG tags")
parser.add_argument("-outdir", help="Tthe output directory for the splitted bams, default is current dir")
parser.add_argument("-barcodes", help='a gz compressed file with accepted barcodes. Optional; if not provided, all cells will be retrieved.')

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

# read only accepted barcodes barcodes
if len(args.barcodes) > 0:
    keep = set()
    with gzip.open(args.barcodes, 'rt') as f:
        i = 0
        for line in f:
            current = 'cellranger_standard:0:1:cb:' + line.strip()
            keep.add(current)
            i += 1
        if i % 100 == 0:
            print('Read %s cells' % i)

print('Found %s accepted barcodes, processing' % len(keep))

## dict with keys are barcode, values are outbam
fouts_dict = {}

prev_cell_barcode = 'empty'
accepted = True
for read in fin:
    tags = read.tags
    RG_list = [ x for x in tags if x[0] == "RG"]
    if RG_list:
        cell_barcode = RG_list[0][1]
        
        if len(args.barcodes) > 0:
            if cell_barcode in keep:
                accepted = True
            else:
                accepted = False
        else:
            accepted = True

        if accepted and cell_barcode not in fouts_dict:
            print('writing %s bamfile' % cell_barcode)
            fout_name = cell_barcode + ".bam"
            fouts_dict[cell_barcode] = pysam.AlignmentFile(os.path.join(args.outdir,fout_name), "wb", template = fin)
            prev_cell_barcode = cell_barcode
        if accepted and cell_barcode in fouts_dict:
            fouts_dict[cell_barcode].write(read)

        ## assumes records are sorted by barcode, so if the current one is different from
        ## the previous one, better close the filehandle
        if accepted and (prev_cell_barcode is not cell_barcode) and (prev_cell_barcode is not 'empty'):
            fouts_dict[prev_cell_barcode].close()
            # print(fouts_dict[prev_cell_barcode].closed)
            prev_cell_barcode = cell_barcode
        
    else: 
        continue
    
fin.close()
# fout.close()
for fout in fouts_dict.values():
    fout.close()
