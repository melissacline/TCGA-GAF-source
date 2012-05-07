#!/usr/bin/env python2.7

import argparse
from BCBio import GFF
import re

parser = argparse.ArgumentParser()
parser.add_argument('inputGff', type=str,  help="Input GFF file")
parser.add_argument('-t', dest="type", default="^.*$", help="Type to output")
args = parser.parse_args()

gffIter = GFF.parse(args.inputGff)
for chrom in gffIter:
    chromName = "chr" + chrom.id
    for hit in chrom.features:
        chromStart = hit.location.nofuzzy_start
        chromEnd = hit.location.nofuzzy_end
        if hit.strand:
            strand = "+"
        else:
            strand = "-"
        print "%s\t%s\t%s\t%s\t1\t%s\t%s\t%s\t0\t1\t%s\t0" \
        % (chromName, chromStart, chromEnd,  hit.qualifiers['Name'][0], strand,
           chromStart, chromEnd, str(int(chromEnd) - int(chromStart) + 1))
