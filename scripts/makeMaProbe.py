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
parser.add_argument('grch37LiteBed', type=str, 
                    help="Input bed file")
parser.add_argument('hg19Bed', type=str, 
                    help="Input bed file")
parser.add_argument("-n", dest="entryNumber", help="Initial entry number",
                    default=0)
args = parser.parse_args()

#
# First, read the hg19 coordiantes
hg19Coords = dict()
hg19CoordsFp = open(args.hg19Bed, "rU")
for line in hg19CoordsFp:
    hg19Bed = Bed.Bed(line.rstrip().split())
    hg19Coords[hg19Bed.name] = hg19Bed
hg19CoordsFp.close()


entryNumber = args.entryNumber
fp = open(args.grch37LiteBed)
for line in fp:
    line = line.rstrip()
    grch37LiteBed = Bed.Bed(line.split())
    hg19Bed = hg19Coords[grch37LiteBed.name]
    gg = Gaf.GafMaProbe(grch37LiteBed, entryNumber)
    geneXrefQuery = """SELECT geneName, locus FROM gafGeneXref
                         WHERE chrom = '%s' AND chromStart <= %s
                           AND chromEnd >= %s
                           and strand = '%s'""" % (hg19Bed.chrom, hg19Bed.chromEnd,
                                                   hg19Bed.chromStart,
                                                   hg19Bed.strand)
    cursor.execute(geneXrefQuery)
    if cursor.rowcount == 1:
        row = cursor.fetchone()
        gg.gene = row["geneName"]
        gg.geneLocus = row["locus"]
    entryNumber = entryNumber + 1
    gg.entryNumber = entryNumber
    gg.write(sys.stdout)
exit(entryNumber)
