#!/usr/bin/env python

import argparse
from pycbio.hgdata import Bed
import Gaf
import MySQLdb
import MySQLdb.cursors
import sys

parser = argparse.ArgumentParser()
parser.add_argument('inputBed', type=str, 
                    help="Input bed file")
parser.add_argument("-n", dest="entryNumber", help="Initial entry number",
                    default=0)
parser.add_argument("-t", dest="exonType", default="",
                    help="Exon type (composite or component)")
args = parser.parse_args()

db = MySQLdb.connect(host="localhost", db="hg19", user="hgcat",
                     passwd="S3attl3-S7u")
cursor = db.cursor(MySQLdb.cursors.DictCursor)

entryNumber = args.entryNumber
fp = open(args.inputBed)
for line in fp:
    line = line.rstrip()
    bb = Bed.Bed(line.split())
    tokens = bb.name.split(";")
    clusterId = tokens.pop()
    bb.name = ";".join(tokens)
    gg = Gaf.GafExon(bb, entryNumber, args.exonType)
    geneXrefQuery = """SELECT geneName, grch37LiteLocus FROM gafGeneXref
                        WHERE clusterId = '%s'""" % (clusterId)
    cursor.execute(geneXrefQuery)
    if cursor.rowcount == 1:
        row = cursor.fetchone()
        gg.gene = row["geneName"]
        gg.geneGrch37LiteLocus = row["grch37LiteLocus"]
        entryNumber = entryNumber + 1
        gg.write(sys.stdout)
exit(entryNumber)
