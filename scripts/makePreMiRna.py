#!/usr/bin/env python2.7

import argparse
import Bio.SeqIO
from pycbio.hgdata import Bed
import Gaf
import MySQLdb
import MySQLdb.cursors
import re
import sys

def getGeneName(xRefList, cursor):
    """Given a list of cross-references, try to find the official HUGO symbol & Entrez Gene ID"""
    #
    # First, try to find an HGNC entry in the xref list
    for xx in xRefList:
        if re.match("^HGNC:", xx):
            cursor.execute("SELECT symbol, entrezId FROM hgnc where hgncId = '%s'" % (xx))
            assert cursor.rowcount <= 1
            if cursor.rowcount == 1:
                row = cursor.fetchone()
                geneName = "%s|%s" % (row['symbol'], row['entrezId'])
                return(geneName)
    #
    # Failing that, try to find an ENTREZGENE entry and search on that.
    for xx in xRefList:
        if re.match("^ENTREZGENE:", xx):
            entrezId = xx.split(":")[1]
            cursor.execute("SELECT symbol FROM hgnc where entrezId = %s" % (entrezId))
            assert cursor.rowcount <= 1
            if cursor.rowcount == 1:
                row = cursor.fetchone()
                geneName = "%s|%s" % (row['symbol'], entrezId)
                return(geneName)
    return("")

            
        
db = MySQLdb.connect(host="localhost", db="hg19", user="hgcat",
                     passwd="S3attl3-S7u")
cursor = db.cursor(MySQLdb.cursors.DictCursor)

parser = argparse.ArgumentParser()
parser.add_argument('preMiRnaBed', type=str, help="Input bed file")
parser.add_argument('miRnaDat', type=str, help="miRNA.dat file from miRBase")
parser.add_argument("-n", dest="entryNumber", help="Initial entry number",
                    default=0)
args = parser.parse_args()
entryNumber = args.entryNumber


#
# First, read the annotation data from miRNA.dat
miRnaDat = dict()
miRnaDatHandle = open(args.miRnaDat, "rU")
for record in Bio.SeqIO.parse(miRnaDatHandle, "embl"):
    miRnaDat[record.name] = record
miRnaDatHandle.close()

#
# Next, read the bed file containing the GRCh37-lite coordinates.
# Convert each line to GAF format and add annotation data from
# miRna.dat
preMiRnaBedFp = open(args.preMiRnaBed)
for line in preMiRnaBedFp:
    bb = Bed.Bed(line.rstrip().split())
    gg = Gaf.GafPreMiRna(bb)
    assert miRnaDat.has_key(bb.name)
    mm = miRnaDat[bb.name]
    gg.featureId = "%s|%s" % (bb.name, mm.id)
    gg.gene = getGeneName(mm.dbxrefs, cursor)
    gg.geneLocus = gg.compositeCoordsToLocus()
    entryNumber = entryNumber + 1
    gg.entryNumber = entryNumber
    gg.write(sys.stdout)
exit(0)
