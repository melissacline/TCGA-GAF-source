#!/usr/bin/env python

import argparse
import Gaf
import re
import sys

def parseFeatureInfo(featureInfo):
    """Given a dbSNP featureInfo record, parse out allele, class, function"""
    (alleleInfo, classInfo, functionInfo) = featureInfo.split(";")
    assert re.search("^Alleles=", alleleInfo)
    assert re.search("^dbSNPclass=", classInfo)
    assert re.search("^dbSNPfunction=", functionInfo)
    alleleData = re.split("^Alleles=", alleleInfo)[1]
    classData = re.split("^dbSNPclass=", classInfo)[1]
    functionData = re.split("^dbSNPfunction=", functionInfo)[1]
    return((alleleData, classData, functionData))

def mergeFeatureInfo(prevFeatureInfo, nextFeatureInfo):
    """Given two GAF dbSNP entries, merge their featureInfo records"""
    (prevAlleles, prevClass, prevFunction) = parseFeatureInfo(prevFeatureInfo)
    (nextAlleles, nextClass, nextFunction) = parseFeatureInfo(nextFeatureInfo)
    newFeatureInfo = "Alleles=%s,%s;dbSNPclass=%s,%s;dbSNPfunction=%s,%s" \
                     % (prevAlleles, nextAlleles, prevClass, nextClass,
                        prevFunction, nextFunction)
    return(newFeatureInfo)


parser = argparse.ArgumentParser()
args = parser.parse_args()

prevGaf = Gaf.Gaf()
prevSnpId = ""
for line in sys.stdin.readlines():
    nextGaf = Gaf.Gaf()
    nextGaf.setFields(line.rstrip().split("\t"))
    if nextGaf.featureId != prevSnpId:
        if len(prevGaf.featureId) > 0:
            prevGaf.write(sys.stdout)
        prevGaf = nextGaf
        prevSnpId = nextGaf.featureId
    else:
        prevGaf.featureCoordinates = "%s;%s" % (prevGaf.featureCoordinates, 
                                                nextGaf.featureCoordinates)
        prevGaf.compositeCoordinates = "%s;%s" % (prevGaf.compositeCoordinates,
                                                  nextGaf.compositeCoordinates)
        prevGaf.gene = "%s;%s" % (prevGaf.gene, nextGaf.gene)
        prevGaf.geneLocus = "%s;%s" % (prevGaf.geneLocus, nextGaf.geneLocus)
        prevGaf.featureInfo = mergeFeatureInfo(prevGaf.featureInfo,
                                               nextGaf.featureInfo)
if nextGaf.featureId != prevGaf.featureId:
    nextGaf.write(sys.stdout)
