#!/usr/bin/env python

import argparse
import Gaf
import sys

    
def lineToGaf(line):
    """Given a line from a gaf file, parse it into a gaf object"""
    tokens = line.rstrip().split('\t')
    gg = Gaf.Gaf()
    gg.setFields(tokens)
    return(gg)


parser = argparse.ArgumentParser()
parser.add_argument('featureGaf', type=str, 
                    help="GAF containing features mapped against the genome")
parser.add_argument('compositeGaf', type=str, 
                    help="GAF containing composites mapped against the genome")
parser.add_argument("-n", dest="entryNumber", help="Initial entry number",
                    default=0)
args = parser.parse_args()

entryNumber = args.entryNumber
#
# Read through the composites.  Make a dictionary with a list of composites
# for each gene.
geneToComposite = dict()
compositeFp = open(args.compositeGaf)
for line in compositeFp:
    comp = lineToGaf(line)
    if geneToComposite.has_key(comp.gene):
        geneToComposite[comp.gene].append(comp)
    else:
        geneToComposite[comp.gene] = [comp]
compositeFp.close()

featureFp = open(args.featureGaf)
for line in featureFp:
    feat = lineToGaf(line)
    compositesThisFeature = geneToComposite[feat.gene]
    for comp in compositesThisFeature:
        featToComp = Gaf.FeatureToCompositeGaf()
        featToComp.assign(feat, comp)
        if len(featToComp.featureCoordinates) > 0:
            featToComp.entryNumber = entryNumber
            print featToComp
            entryNumber = entryNumber + 1
