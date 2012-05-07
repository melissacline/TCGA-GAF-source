#!/usr/bin/env python2.7

import argparse
import Gaf
import sys



prevGaf = Gaf.Gaf()
prevFeatureId = ""
for line in sys.stdin.readlines():
    nextGaf = Gaf.Gaf()
    nextGaf.setFields(line.rstrip().split("\t"))
    if nextGaf.featureId != prevFeatureId:
        if len(prevGaf.featureId) > 0:
            prevGaf.write(sys.stdout)
        prevGaf = nextGaf
        prevFeatureId = nextGaf.featureId
    else:
        prevGaf.featureCoordinates = "%s;%s" % (prevGaf.featureCoordinates, 
                                                nextGaf.featureCoordinates)
        prevGaf.compositeCoordinates = "%s;%s" % (prevGaf.compositeCoordinates,
                                                  nextGaf.compositeCoordinates)
        prevGaf.gene = "%s;%s" % (prevGaf.gene, nextGaf.gene)
        prevGaf.geneLocus = "%s;%s" % (prevGaf.geneLocus, nextGaf.geneLocus)
        prevGaf.featureInfo = "%s;%s" % (prevGaf.featureInfo, nextGaf.featureInfo)
if nextGaf.featureId != prevGaf.featureId:
    nextGaf.write(sys.stdout)
