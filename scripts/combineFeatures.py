#!/usr/bin/env python

import argparse
import Gaf
import re
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
        if not re.search(re.sub("\?", "\?", nextGaf.gene), prevGaf.gene):
            prevGaf.gene = "%s;%s" % (prevGaf.gene, nextGaf.gene)
        if not re.search(re.sub("\+", "\+", nextGaf.geneLocus), prevGaf.geneLocus):
            prevGaf.geneLocus = "%s;%s" % (prevGaf.geneLocus, nextGaf.geneLocus)
        if not re.search(nextGaf.featureInfo, prevGaf.featureInfo):
            prevGaf.featureInfo = "%s;%s" % (prevGaf.featureInfo, nextGaf.featureInfo)
if nextGaf.featureId == prevGaf.featureId:
    prevGaf.write(sys.stdout)

