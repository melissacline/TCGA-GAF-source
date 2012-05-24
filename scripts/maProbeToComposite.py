#!/usr/bin/env python

import argparse
from pycbio.hgdata import Bed
import Gaf
import MySQLdb
import MySQLdb.cursors
import re
import sys

def findOverlappingMaProbes(compositeBed, cursor):
    """Given bed coordinates for a composite, return the names of any maProbe
    that overlaps it. If no maProbe overlaps it, return an empty list"""
    overlappingProbeList = list()
    query = """SELECT name FROM gafMaProbe WHERE chrom = '%s'
                  AND chromStart <= %s AND chromEnd >= %s
                  AND strand = '%s'""" % (compositeBed.chrom,
                                          compositeBed.chromEnd,
                                          compositeBed.chromStart,
                                          compositeBed.strand)
    cursor.execute(query)
    print query
    for row in cursor.fetchall():
        overlappingProbeList.append(row['name'])
    return overlappingProbeList

db = MySQLdb.connect(host="localhost", db="hg19", user="hgcat",
                     passwd="S3attl3-S7u")
cursor = db.cursor(MySQLdb.cursors.DictCursor)
            
parser = argparse.ArgumentParser()
parser.add_argument('maProbeGaf', type=str, help="maProbe GAF file")
parser.add_argument('compositeGaf', type=str, help="composite GAF file")
parser.add_argument('compositeBed', type=str,
                    help="composite BED file, hg19 coordinates")
parser.add_argument("-n", dest="entryNumber", help="Initial entry number",
                    default=0)
args = parser.parse_args()

entryNumber = args.entryNumber

#
# First, read the maProbe GAF data into a dictionary.  
maProbe = dict()
maProbeGafFp = open(args.maProbeGaf)
for line in maProbeGafFp:
    maProbeGaf = Gaf.Gaf()
    maProbeGaf.setFields(line.rstrip().split("\t"))
    maProbe[maProbeGaf.featureId] = maProbeGaf
maProbeGafFp.close()

#
# Next, read the composite bed data into a dictionary.  
compositeBed = dict()
compositeBedFp = open(args.compositeBed)
for line in compositeBedFp:
    bedThisComposite = Bed.Bed(line.rstrip().split())
    bedName = bedThisComposite.name.split("|")[0]
    compositeBed[bedName] = bedThisComposite
compositeBedFp.close()

#
# Next, read the composite GAF file.  For each composite, use the
# hg19 bed coordinates to find any overlapping maProbe records.
# If one is found, look up the maProbe gaf from the dictionary,
# and output the mapping.
compositeGafFp = open(args.compositeGaf)
for line in compositeGafFp:
    compositeGaf = Gaf.Gaf()
    compositeGaf.setFields(line.rstrip().split("\t"))
    featureId = compositeGaf.featureId.split("|")[0]
    assert compositeBed.has_key(featureId)
    compositeBedEntry = compositeBed[featureId]
    overlappingMaProbeNames = findOverlappingMaProbes(compositeBedEntry, cursor)
    for maProbeName in overlappingMaProbeNames:
        assert maProbe.has_key(maProbeName)
        maProbeGaf = maProbe[maProbeName]
        maProbeToCompositeGaf = Gaf.FeatureToCompositeGaf()
        maProbeToCompositeGaf.assign(maProbeGaf, compositeGaf)
        if len(maProbeToCompositeGaf.featureCoordinates) > 0:
            entryNumber = entryNumber + 1
            maProbeToCompositeGaf.entryNumber = entryNumber
            maProbeToCompositeGaf.write(sys.stdout)
exit(0)
