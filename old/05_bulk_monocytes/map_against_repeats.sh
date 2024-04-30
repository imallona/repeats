#!/bin/bash
#
# Izaskun Mallona
#
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80095 against human repeats
#
# 24 aug 2020
# requires get_data.sh to be run first

# maybe snakemake made?

snakemake -j 20 -s bulk_rnaseq.snmk
