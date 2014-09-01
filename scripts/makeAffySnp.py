#!/usr/bin/env python

import argparse
from pycbio.hgdata import Bed
import Gaf
import Grch37LiteGaf
import MySQLdb
import MySQLdb.cursors
import sys

parser = argparse.ArgumentParser()
parser.add_argument('affySnpBed', type=str, 
                    help="Affy SNP data in GRCh37-lite coordinates")
parser.add_argument('oldAffySnpGaf', type=str, 
                    help="Gaf 2.1 file with AffySnp-Genome records")
parser.add_argument("-n", dest="entryNumber", help="Initial entry number",
                    default=0)
args = parser.parse_args()

db = MySQLdb.connect(host="localhost", db="hg19", read_default_file="~/.my.cnf")
cursor = db.cursor(MySQLdb.cursors.DictCursor)

entryNumber = int(args.entryNumber)

#
# Make a dictionary of the old entries.
oldAffySnp = dict()
oldAffySnpGafFp = open(args.oldAffySnpGaf)
for line in oldAffySnpGafFp:
    oldAffySnpGaf = Gaf.Gaf(line)
    oldAffySnp[oldAffySnpGaf.featureId] = oldAffySnpGaf
oldAffySnpGafFp.close()

fp = open(args.affySnpBed)
for line in fp:

    #
    # Get the basic fields from the BED data
    bb = Bed.Bed(line.split())
    gg = Grch37LiteGaf.GafAffySnp(bb, entryNumber)

    #
    # From the old GAF data, get the featureInfo string
    gg.featureInfo = oldAffySnp[gg.featureId].featureInfo

    #
    # Look up the gene and gene locus for these coordinates
    gg.gene = ""
    gg.geneLocus = ""
    geneXrefQuery = """SELECT geneName, grch37LiteLocus FROM gafGeneXref
                        WHERE grch37LiteChrom = '%s'
                          AND grch37LiteChromStart <= %d
                          AND grch37LiteChromEnd >= %d""" % (bb.chrom,
                                                             bb.chromEnd,
                                                             bb.chromStart)
    cursor.execute(geneXrefQuery)
    delimiter=""
    for row in cursor.fetchall():
        gg.gene = "%s%s%s" % (gg.gene, delimiter, row["geneName"])
        gg.geneLocus = "%s%s%s" % (gg.geneLocus, delimiter,
                                   row["grch37LiteLocus"])
        delimiter = ","
    entryNumber = entryNumber + 1
    gg.write(sys.stdout)
exit(0)
