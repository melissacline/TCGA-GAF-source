#!/usr/bin/env python

import argparse
from pycbio.hgdata import Bed
import Gaf
import MySQLdb
import MySQLdb.cursors
import sys

db = MySQLdb.connect(host="localhost", db="hg19", user="hgcat",
                     passwd="S3attl3-S7u")
cursor = db.cursor(MySQLdb.cursors.DictCursor)

parser = argparse.ArgumentParser()
parser.add_argument('inputBed', type=str, 
                    help="Input bed file")
parser.add_argument("-n", dest="entryNumber", help="Initial entry number",
                    default=0)
args = parser.parse_args()

entryNumber = args.entryNumber
fp = open(args.inputBed)
for line in fp:
    line = line.rstrip()
    bb = Bed.Bed(line.split())
    tokens = bb.name.split(";")
    clusterId = tokens.pop()
    bb.name = ";".join(tokens)
    gg = Gaf.GafTranscript(bb, entryNumber)
    kgXrefQuery = "SELECT value FROM knownToRefSeq WHERE name = '%s'" \
                  % (bb.name)
    cursor.execute(kgXrefQuery)
    if cursor.rowcount == 1:
        row = cursor.fetchone()
        if len(row["value"]) > 0:
            gg.featureAliases = row["value"]
    geneXrefQuery = """SELECT geneName, locus FROM gafGeneXref
                        WHERE clusterId = '%s'""" % (clusterId)
    cursor.execute(geneXrefQuery)
    if cursor.rowcount == 1:
        row = cursor.fetchone()
        gg.gene = row["geneName"]
        gg.geneLocus = row["locus"]
        entryNumber = entryNumber + 1
        gg.write(sys.stdout)
exit(entryNumber)
