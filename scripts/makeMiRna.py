#!/usr/bin/env python2.7

import argparse
import Bio.SeqIO
from pycbio.hgdata import Bed
import Gaf
import re
import sys

            

parser = argparse.ArgumentParser()
parser.add_argument('preMiRnaGaf', type=str, help="pre-miRNA GAF file",
                    default="data/gaf/preMiRna.gaf")
parser.add_argument('inputDat', type=str, help="miRNA.dat file from miRBase",
                    default="data/miRna.dat")
parser.add_argument("-n", dest="entryNumber", help="Initial entry number",
                    default=0)
args = parser.parse_args()

entryNumber = args.entryNumber

#
# First, read the annotation data from miRNA.dat
miRnaDat = dict()
miRnaDatHandle = open(args.inputDat, "rU")
for record in Bio.SeqIO.parse(miRnaDatHandle, "embl"):
    miRnaDat[record.name] = record
miRnaDatHandle.close()


#
# Next, read the pre-miRNA GAF data.  For each pre-miRNA,
# check the features data from the miRnaDat dictionary
# to identify any miRnas.  Compute the genomic coordinates
# of the miRnas, and output with appropriate annotations.
preMiRnaFp = open(args.preMiRnaGaf)
for line in preMiRnaFp:
    preMiRnaGaf = Gaf.Gaf()
    preMiRnaGaf.setFields(line.rstrip().split("\t"))
    preMiRnaName = preMiRnaGaf.featureId.split("|")[0]
    assert miRnaDat.has_key(preMiRnaName)
    for miRnaFeature in miRnaDat[preMiRnaName].features:
        if miRnaFeature.type == 'miRNA':
            miRnaStart = miRnaFeature.location.nofuzzy_start + 1
            miRnaEnd = miRnaFeature.location.nofuzzy_end
            if miRnaFeature.qualifiers.has_key('accession'):
                accession = miRnaFeature.qualifiers['accession'][0]
                product = miRnaFeature.qualifiers['product'][0]
            else:
                accession = ""
                product = ""
            miRnaGaf = Gaf.Gaf()
            miRnaGaf.copy(preMiRnaGaf)
            miRnaGaf.featureId = "%s|%s" % (product, accession)
            miRnaGaf.featureType = "miRNA"
            miRnaGaf.compositeCoordinates \
                = miRnaGaf.subsetCompositeCoordinates(miRnaStart, miRnaEnd)
            miRnaGaf.featureCoordinates = "1-%d" % (miRnaEnd - miRnaStart + 1)
            miRnaGaf.featureInfo = "pre-miRNA=%s" % (preMiRnaGaf.featureId)
            entryNumber = entryNumber + 1
            miRnaGaf.entryNumber = entryNumber
            miRnaGaf.write(sys.stdout)
exit(0)
