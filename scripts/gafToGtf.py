#!/usr/bin/env python

import sys, argparse
from GafParse import Transcript

parser = argparse.ArgumentParser(description="Create GTF file from input GAF file")
parser.add_argument('inputGaf', type=str,help="Input GAF file")
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()

f = open(args.inputGaf,'r')

for line in f:
    tx = Transcript(line)
    if tx.FeatureType == 'junction':
        tx.makeJunctionGTF()
    else:
        tx.makeExonGTF()
    for l in tx.exonGTF:
        print l
    if tx.FeatureType == 'transcript' and tx.CompositeType == 'genome':
        tx.mapCDS()
        tx.makeORFGTF()
        for l in tx.orfGTF:
            print l
        print
f.close()

