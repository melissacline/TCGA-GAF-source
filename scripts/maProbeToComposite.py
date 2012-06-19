#!/usr/bin/env python

import argparse
from pycbio.hgdata import Bed
import Gaf
import MySQLdb
import MySQLdb.cursors
import re
import sys

def findOverlappingMaProbes(chrom, chromStart, chromEnd, strand, cursor, debug):
    """Given bed coordinates for a composite, return a list of any maProbe
    that overlaps it. If some probe overlaps it, return a key to match the
    dictionary created below.  If no maProbe overlaps it, return an empty list.
    Note that BED data (the table) is 0-based while GAF data (everything else)
    is 1-based.  Adjust appropriately"""
    overlappingProbeList = list()
    query = """SELECT name, chrom, chromStart, chromEnd, strand FROM gafMaProbeGrch37Lite
                WHERE chrom = '%s' AND chromStart <= %s
                  AND chromEnd >= %d AND strand = '%s'"""  % (chrom, chromEnd,
                                                              int(chromStart) - 1, strand)
    if debug:
        print "executing", query
    cursor.execute(query)
    for row in cursor.fetchall():
        newEntry = "%s:%s:%d-%s:%s" % (row['name'], row['chrom'],
                                       int(row['chromStart']) + 1,
                                       row['chromEnd'], row['strand'])
        if debug:
            print "adding", newEntry, "to results"
        overlappingProbeList.append(newEntry)
    return overlappingProbeList

db = MySQLdb.connect(host="localhost", db="hg19", user="hgcat",
                     passwd="S3attl3-S7u")
cursor = db.cursor(MySQLdb.cursors.DictCursor)
            
parser = argparse.ArgumentParser()
parser.add_argument('maProbeGaf', type=str, help="maProbe GAF file")
parser.add_argument('compositeGaf', type=str, help="composite GAF file")
parser.add_argument("-n", dest="entryNumber", help="Initial entry number",
                    default=0)
parser.add_argument("-d", dest="debug", help="display debugging messages",
                    default=False)
args = parser.parse_args()

entryNumber = args.entryNumber


#
# First, read the maProbe GAF data into a dictionary.  One probe can have
# many alignments, so use a key based on the name and coordinates.
maProbe = dict()
maProbeGafFp = open(args.maProbeGaf)
for line in maProbeGafFp:
    maProbeGaf = Gaf.Gaf(line)
    (maChrom, maCoordinateString, maStrand) = maProbeGaf.compositeCoordinates.split(":")
    maChromStart = maCoordinateString.split("-")[0]
    maChromEnd = maCoordinateString.split("-")[-1]
    key = "%s:%s:%s-%s:%s" % (maProbeGaf.featureId, maChrom, maChromStart, maChromEnd,
                              maStrand)
    maProbe[key] = maProbeGaf
maProbeGafFp.close()


#
# Next, read the composite GAF file.  For each composite, use the
# hg19 bed coordinates to find any overlapping maProbe records.
# If one is found, look up the maProbe gaf from the dictionary,
# and output the mapping.
compositeGafFp = open(args.compositeGaf)
for line in compositeGafFp:
    compositeGaf = Gaf.Gaf(line)
    if args.debug:
        print "working on", compositeGaf
    (cgChrom, cgCoordinateString, cgStrand) = compositeGaf.compositeCoordinates.split(":")
    cgCoordStart = cgCoordinateString.split("-")[0]
    cgCoordEnd = cgCoordinateString.split("-")[-1]
    overlappingMaProbeNames = findOverlappingMaProbes(cgChrom, cgCoordStart, cgCoordEnd,
                                                      cgStrand, cursor, args.debug)
    for maProbeEntry in overlappingMaProbeNames:
        if args.debug:
            print "workng on probe", maProbeEntry
        maProbeGaf = maProbe[maProbeEntry]
        maProbeToCompositeGaf = Gaf.FeatureToCompositeGaf()
        maProbeToCompositeGaf.assign(maProbeGaf, compositeGaf, mergeAdjacent=True,
                                     debug=args.debug)
        if len(maProbeToCompositeGaf.featureCoordinates) > 0:
            #
            # Unlike other feature-composite mappings, this time we'll take the
            # gene and gene locus fields from the composite instead of the
            # feature (i.e. instead of the probe)
            maProbeToCompositeGaf.gene = compositeGaf.gene
            maProbeToCompositeGaf.geneLocus = compositeGaf.geneLocus
            entryNumber = entryNumber + 1
            maProbeToCompositeGaf.entryNumber = entryNumber
            maProbeToCompositeGaf.write(sys.stdout)
exit(0)
