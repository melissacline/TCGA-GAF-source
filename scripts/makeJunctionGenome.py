#!/usr/bin/env python

import argparse
from pycbio.hgdata import Bed
import Grch37LiteGaf
import MySQLdb
import MySQLdb.cursors
import sys

db = MySQLdb.connect(host="localhost", db="hg19", read_default_file="~/.my.cnf")
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
    geneXrefQuery = """SELECT DISTINCT geneName, grch37LiteLocus
                         FROM gafGeneXref
                        WHERE clusterId = '%s'""" % (clusterId)
    cursor.execute(geneXrefQuery)
    if cursor.rowcount == 1:
        row = cursor.fetchone()
        gene = row["geneName"]
        geneLocus = row["grch37LiteLocus"]
        for ii in range(len(bb.blocks)-1):
            gg = Grch37LiteGaf.GafJunction(bb, entryNumber, ii)
            gg.gene = gene
            gg.geneLocus = geneLocus
            gg.write(sys.stdout)
            entryNumber = entryNumber + 1
exit(entryNumber)
