#!/usr/bin/env python2.7

import argparse
from pycbio.hgdata import Bed
import Gaf
import sys

parser = argparse.ArgumentParser()
parser.add_argument('inputBed', type=str, 
                    help="Input bed file")
parser.add_argument("-n", dest="entryNumber", help="Initial entry number",
                    default=0)
args = parser.parse_args()

entryNumber = args.entryNumber
fp = open(args.inputBed)
for line in fp:
    line = line.rstrip()
    bb = Bed.Bed(line.split())
    gg = Gaf.GafAffySnp(bb, entryNumber)
    entryNumber = entryNumber + 1
    gg.write(sys.stdout)
exit(entryNumber)
