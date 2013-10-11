#!/usr/bin/env python

import argparse
import Gaf
import re

def rangeCleanup(blockList):
    """Given a set of coordinate range blocks, put them into proper order
       w.r.t. the plus strand, such that start < end for each block and
       block[i].start < block[i+1].start for all i"""
    # Run through the list and reorder the start and end coordinates
    # of each block if needed.
    for ii in range(len(blockList)):
        if re.search("-", blockList[ii]):
            start = int(blockList[ii].split("-")[0])
            end = int(blockList[ii].split("-")[1])
            if start > end:
                blockList[ii] = "%d-%d" % (end, start)
    # Now run through the list and sort the blocks into increasing order
    listSorted = False
    while not listSorted:
        swappedPositions = False
        for ii in range(len(blockList) - 1):
            start1 = int(blockList[ii].split("-")[0])
            start2 = int(blockList[ii+1].split("-")[0])
            if start1 > start2:
                temp = blockList[ii]
                blockList[ii] = blockList[ii+1]
                blockList[ii+1] = temp
                swappedPositions = True
        if not swappedPositions:
            listSorted = True
    return(blockList)
        

def gafLineToBedLine(gafData):
    """Given a gaf entry, extract the bed coordinates"""
    allBedLines = "";
    bedLineDelimiter = "";
    name = gafData.featureId
    for subCoordinates in gafData.compositeCoordinates.split(";"):
        if re.search(":", subCoordinates):
            (chrom, coordinateBlocks, strand) = subCoordinates.split(":")
            chrom = re.sub("chrM_rCRS", "chrM", chrom)
            coordinateRanges = rangeCleanup(coordinateBlocks.split(","))
            if re.search("-", coordinateRanges[0]):
                start = int(coordinateRanges[0].split("-")[0]) - 1
                end = int(coordinateRanges[-1].split("-")[-1])
            else:
                start = int(coordinateRanges[0]) - 1
                end = int(coordinateRanges[0])
            if strand == '-' and start > end:
                temp = start + 1
                start = end - 1
                end = temp
            blockStarts = ""
            blockLengths = ""
            delimiter = ""
            for ii in range(len(coordinateRanges)):
                if re.search("-", coordinateRanges[ii]):
                    thisStart = int(coordinateRanges[ii].split("-")[0]) - 1
                    thisEnd = int(coordinateRanges[ii].split("-")[-1])
                else:
                    thisStart = int(coordinateRanges[ii]) - 1
                    thisEnd = int(coordinateRanges[ii])
                if strand == "-" and thisStart > thisEnd:
                    temp = thisStart + 1
                    thisStart = thisEnd - 1
                    thisEnd = temp
                blockLengths = blockLengths + delimiter \
                               + str(thisEnd - thisStart)
                blockStarts = blockStarts + delimiter + str(thisStart - start)
                delimiter = ","
            bedLine = "\t".join((chrom, str(start), str(end), name, "1",
                                 strand, str(start), str(end), "0",
                                 str(len(coordinateRanges)), blockLengths,
                                 blockStarts))
            allBedLines = "%s%s%s" % (allBedLines, bedLineDelimiter, bedLine)
            bedLineDelimiter = "\n"
    return(allBedLines)
                       

    


parser = argparse.ArgumentParser("gaf to bed")
parser.add_argument('inputGaf', type=str, help="Input GAF file")
parser.add_argument("-d", dest="debug", help="Turn on debugging", default=False)
args = parser.parse_args()

fp = open(args.inputGaf)
for row in fp:
    row = row.rstrip()
    if args.debug:
        print "Input Row", row
    gafData = Gaf.Gaf(row)
    bedLine = gafLineToBedLine(gafData)
    if bedLine != "":
        print bedLine
    if args.debug:
        print "---"
