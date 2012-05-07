#!/usr/bin/env python2.7

import argparse
from pycbio.hgdata import Bed
from pycbio.hgdata import Psl
import re

def bedIsSubsetOfPsl(bedAlign, pslAlign):
    """Return True if the bed alignment falls within the bounds of the PSL align"""
    if bedAlign.chrom == pslAlign.tName and bedAlign.strand == pslAlign.strand:
        if pslAlign.tStart <= bedAlign.chromEnd \
               and pslAlign.tEnd >= bedAlign.chromStart:
            return(True)
    return(False)


parser = argparse.ArgumentParser()
parser.add_argument('inputPsl', type=str,
                    help="PSL file from which selected records will be printed") 
parser.add_argument('inputBed', type=str, 
                    help="Bed file to be used for screening the PSL records")
args = parser.parse_args()

bedAlignments = dict()
bedFp = open(args.inputBed)
for line in bedFp:
    line = line.rstrip()
    bedAlign = Bed.Bed(line.split())
    bedAlignments[bedAlign.name] = bedAlign
bedFp.close()
pslFp = open(args.inputPsl)
readPastHeader = False
for line in pslFp:
    line = line.rstrip()
    if re.search("^--", line):
        readPastHeader = True
    else:
        if readPastHeader:
            pslAlign = Psl.Psl(line.split())
            if bedAlignments.has_key(pslAlign.qName):
                if bedIsSubsetOfPsl(bedAlignments[pslAlign.qName], pslAlign):
                    print line
