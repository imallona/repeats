#!/bin/python
#
# https://github.com/pezmaster31/bamtools/issues/135
# jvhaarst commented on Feb 16, 2018

import os
import fileinput
for line in fileinput.input():
    line=line.strip()
    tags = line.split()[11:]
    for tag in tags:
        tag_split=tag.split(':')
        if 'CB' in tag_split:
            barcode=tag_split[-1]
            directory = barcode[:5]
            outfile = directory+"/"+barcode
            try:
                f=open(outfile,"a+")
                print(line, file=f)
                f.close()
            except:
                os.mkdir(directory)
                f=open(outfile,"a+")
                print(line, file=f)
                f.close()
