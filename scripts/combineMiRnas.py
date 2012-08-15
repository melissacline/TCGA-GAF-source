#!/usr/bin/env python
"""
combineMiRnas.py: combine the GAf entries for miRNAs with the same feature ID

Usage:
combimeMiRnas.py uncombined.gaf > combined.gaf

Given a GAF file, in which each line represents a single miRNA genomic alignment,
generate a GAF file in which each distinct miRNA has one line containing all of
its genomic alignments, with feature and composite coordinate strings delimited
by semicolons and featureInfo (pre-miRNA=) entries combined and delimimted with
commas.

Assumptions:
The input data should be sorted by miRNA feature ID
"""

import argparse
import Grch37LiteGaf
import re
import sys

parser = argparse.ArgumentParser()
parser.add_argument('inputGaf', type=str, 
                    help="Input GAF file: if any feature has two alignments, combine the records")
args = parser.parse_args()

prevRnaId = ""
prevAccession = ""
fp = open(args.inputGaf)
for line in fp:
    nextGaf = Grch37LiteGaf.GafMiRna(line, createFromBedInput=False)
    (nextRnaId,nextAccession) = nextGaf.featureId.split("|")
    nextRnaIdTokens = nextRnaId.split("-")
    if len(nextRnaIdTokens) == 4:
        nextRnaId = "%s-%s-%s" % (nextRnaIdTokens[0], nextRnaIdTokens[1],
                                  nextRnaIdTokens[2])
    if nextRnaId != prevRnaId or nextAccession != prevAccession:
        if len(prevRnaId) > 0 and len(prevAccession) > 0:
            prevGaf.write(sys.stdout)
        prevGaf = nextGaf
        prevRnaId = nextRnaId
        prevAccession = nextAccession
    else:
        prevGaf.combine(nextGaf)
if nextGaf.featureId == prevGaf.featureId:
    prevGaf.write(sys.stdout)
