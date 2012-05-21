#!/usr/bin/env python

import argparse
import Bio.SeqIO
from pycbio.hgdata import Bed
import Gaf
from BCBio import GFF
import re
import sys

import re

parser = argparse.ArgumentParser()
parser.add_argument('miRnaBed', type=str, help="Input bed file")
parser.add_argument('inputGff', type=str,  help="Input GFF file")
parser.add_argument("-n", dest="entryNumber", help="Initial entry number",
                    default=0)
args = parser.parse_args()

entryNumber = args.entryNumber

#
# Read over the gff file collecting data for ID mapping. 
# map the ID to what will be the displayed ID: the combination of the
# name and accession number.  For miRNAs, also record the ID of the
# pre-miRNA that it was derived from.
miRnaToPreMiRna = dict()
idToLabel = dict()
gffIter = GFF.parse(args.inputGff)
for chrom in gffIter:
    for hit in chrom.features:
        id = hit.id
        label = "%s|%s" % (hit.qualifiers["Name"][0], hit.qualifiers["accession_number"][0])
        idToLabel[id] = label
        if hit.type == "miRNA":
            miRnaToPreMiRna[hit.id] = hit.qualifiers["derives_from"][0]
            

#
# Read the bed file containing the GRCh37-lite coordinates.
# While converting each line to GAF format, replace the ID
# with the miRNA name, look up the name of the pre-miRNA that
# the miRNA is derived from, and note that in the featureInfo field.
miRnaBedFp = open(args.miRnaBed)
for line in miRnaBedFp:
    bb = Bed.Bed(line.rstrip().split())
    gg = Gaf.GafMiRna(bb)
    preMiRnaId = miRnaToPreMiRna[gg.featureId]
    preMiRnaName = idToLabel[preMiRnaId]
    gg.featureInfo = "pre-miRNA=%s" % (preMiRnaName)
    gg.featureId = idToLabel[gg.featureId]
    entryNumber = entryNumber + 1
    gg.entryNumber = entryNumber
    gg.write(sys.stdout)
exit(0)


