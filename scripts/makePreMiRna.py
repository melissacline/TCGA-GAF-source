#!/usr/bin/env python

import argparse
import Bio.SeqIO
from pycbio.hgdata import Bed
import Grch37LiteGaf
import MySQLdb
import MySQLdb.cursors
import re
import sys

def getGeneName(xRefList, id, cursor):
    """Given a list of cross-references, try to find the official HUGO symbol & Entrez Gene ID"""
    #
    # First, try to find an ENTREZGENE entry in the xref list
    for xx in xRefList:
        if re.match("^ENTREZGENE:", xx):
            entrezId = xx.split(":")[1]
            cursor.execute("SELECT symbol FROM hgnc where entrezId = %s" % (entrezId))
            assert cursor.rowcount <= 1
            if cursor.rowcount == 1:
                row = cursor.fetchone()
                geneName = "%s|%s" % (row['symbol'], entrezId)
                return(geneName)
    #
    # Failing that, try to find an HGNC entry and search on that.
    # One might think that the HGNC entry would be the best thing to search on.
    # But it turns out that at this time (May 20, 2012), the Entrez Gene IDs distinguish
    # between paralogs more effectively.
    for xx in xRefList:
        if re.match("^HGNC:", xx):
            cursor.execute("SELECT symbol, entrezId FROM hgnc where hgncId = '%s'" % (xx))
            assert cursor.rowcount <= 1
            if cursor.rowcount == 1:
                row = cursor.fetchone()
                geneName = "%s|%s" % (row['symbol'], row['entrezId'])
                return(geneName)

    #
    # As a last ditch effort, try to find the miRBase ID in the list of synonyms.
    # When doing so, take care to avoid degenerate matches.  Look for an exact
    # match to one of the tokens in the symbol string, in which multiple synonyms
    # are each separated by a comma followed by a single space.  If there is a hit
    # on the synonym, accept it if and only if there is a match to one HGNC entry.
    cursor.execute("""SELECT symbol, entrezId FROM hgnc
                       WHERE synonyms REGEXP '^((.*), )*%s(, (.*))*$'""" % (id))
    if cursor.rowcount == 1:
        row = cursor.fetchone()
        geneName = "%s|%s" % (row['symbol'], row['entrezId'])
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
# First, read the annotation data from miRNA.dat.  For each human
# record (signified by the name beginning with "hsa"), save the
# data in a dictonary indexed by accession 
miRnaDat = dict()
miRnaDatHandle = open(args.miRnaDat, "rU")
for record in Bio.SeqIO.parse(miRnaDatHandle, "embl"):
    if re.search("^hsa", record.name):
        miRnaDat[record.id] = record
miRnaDatHandle.close()


#
# Next, read the bed file containing the GRCh37-lite coordinates.
# Convert each line to GAF format and add annotation data from
# miRna.dat
preMiRnaBedFp = open(args.preMiRnaBed)
for line in preMiRnaBedFp:
    bb = Bed.Bed(line.rstrip().split())
    gg = Grch37LiteGaf.GafPreMiRna(bb)
    assert miRnaDat.has_key(bb.name)
    mm = miRnaDat[bb.name]
    gg.featureId = "%s|%s" % (mm.name, mm.id)
    gg.gene = getGeneName(mm.dbxrefs, id, cursor)
    gg.geneLocus = gg.compositeCoordsToLocus()
    entryNumber = entryNumber + 1
    gg.entryNumber = entryNumber
    gg.write(sys.stdout)
exit(0)
