#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('inputFile1', type=str, help="Input file 1")
parser.add_argument('inputFile2', type=str, help="Input file 2")
parser.add_argument("-d", dest="delimiter", help="Field delimiter",
                    default="\t")
parser.add_argument("-i", dest="ignoreColumns", help="Columns to ignore",
                    default="")
args = parser.parse_args()
columnsToIgnore = set(args.ignoreColumns.split(","))

fp1 = open(args.inputFile1)
fp2 = open(args.inputFile2)
lineCounter = 1
for line1 in fp1:
    line1 = line1.rstrip()
    line2 = fp2.readline().rstrip()
    tokens1 = line1.split(args.delimiter)
    tokens2 = line2.split(args.delimiter)
    allColumnsMatch = True
    if len(tokens1) != len(tokens2):
        print "Line", lineCounter, "token length mismatch: ", len(tokens1), "vs", \
              len(tokens2)
    for ii in range(min(len(tokens1), len(tokens2))):
        if not str(ii+1) in columnsToIgnore:
            if tokens1[ii] != tokens2[ii]:
                print "line", lineCounter, "column", ii+1, "mismatch:", \
                      tokens1[ii], tokens2[ii]
                allColumnsMatch = False
    if len(tokens1) > len(tokens2):
        for ii in range(len(tokens2), len(tokens1)):
            if not str(ii+1) in columnsToIgnore:
                if len(tokens1[ii]) > 0:
                    print "line", lineCounter, "extra column", ii, "in first file:", \
                          tokens1[ii]
                    allColumnsMatch = False
    if len(tokens2) > len(tokens1):
        for ii in range(len(tokens1), len(tokens2)):
            if not str(ii+1) in columnsToIgnore:
                if len(tokens2[ii]) > 0:
                    print "line", lineCounter, "extra column", ii, "in second file:", \
                          tokens2[ii]
                    allColumnsMatch = False
    if allColumnsMatch:
        print "line", lineCounter, "match"
    lineCounter = lineCounter + 1
