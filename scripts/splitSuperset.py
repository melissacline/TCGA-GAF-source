#!/usr/bin/env python

import argparse
import Gaf
import os
import re

parser = argparse.ArgumentParser()
parser.add_argument('supersetGaf', type=str, help="superset GAF file")
args = parser.parse_args()

outputFp = None
prevTargetGafFile = ""
gafFp = open(args.supersetGaf)
for line in gafFp:
    #
    # Parse the GAF entry out of the line, and see what the feature and
    # composite types are.
    gafEntry = Gaf.Gaf(line)
    featureType = gafEntry.featureType
    compositeType = gafEntry.compositeType
    #
    # Assemble the target filename from the feature and composite types.
    # If the target filename is the same as it was for the last name,
    # then print this line into the current file handle.  Otherwise, close
    # the current file handle, open a new file handle to write, and write
    # the line to this file.  This takes advantage of the fact that the GAF
    # file in earlier processing steps is assembled from chunks, and those
    # chunks are defined by feature and composite type.
    targetGafFile = "%s.%s.gaf" % (featureType, compositeType)
    if targetGafFile != prevTargetGafFile:
        if outputFp != None:
            outputFp.close()
        outputFp = open(targetGafFile, "w")
        prevTargetGafFile = targetGafFile
    print >>outputFp, line
outputFp.close()
        
    
