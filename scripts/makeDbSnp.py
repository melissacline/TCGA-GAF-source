#!/usr/bin/env python

import argparse
from pycbio.hgdata import Bed
import Grch37LiteGaf
import MySQLdb
import MySQLdb.cursors
import re
import sys

def basicSnpInfo(bedData, snpTable, cursor):
    """Fill in the basic info from the dbSnp table"""
    dbSnpQuery = """select observed, class, func from %s where name = '%s'
                      and chrom = '%s' and chromStart = %d and chromEnd = %d""" \
                 % (snpTable, bedData.name, bedData.chrom, bedData.chromStart,
                    bedData.chromEnd)
    cursor.execute(dbSnpQuery)
    assert cursor.rowcount == 1
    row = cursor.fetchone()
    snpInfo = "Alleles=%s;dbSNPclass=%s;dbSNPfunction=%s" \
              % (row['observed'], row['class'], row['func'])
    return(snpInfo)

def mapSnpToLocus(bedData, snpInfo, cursor):
    """Determine which locus, if any, overlaps this SNP"""
    #
    # By definition, a SNP is designated as 'near-gene-5' if it is
    # no more than 2000 bases upstream of a gene and/or designated
    # as 'near-gene-3' if it is no more than 500 bases downstream of
    # a gene.  To look for an overlapping gene, prepare to look for something
    # that overlaps the start and end coordinates.  If the designation includes
    # 'near-gene-5', then increase the end coordinate by 2000 bases.  When you look
    # for a gene that overlaps the start-end range, this will effectively look for
    # something at the SNP or up to 2000 bases upstream.  Similarly, if the
    # designation includes 'near-gene-3', then decrease the start coordinate
    # by 500 bases.  Note that the same SNP can be 'near-gene-3' and 'near-gene-5'.
    #
    # Flip all of this around for the minus strand.
    #
    searchStart = bedData.chromStart
    searchEnd = bedData.chromEnd
    if re.search("near-gene-5", snpInfo):
        if bedData.strand == '+':
            searchEnd = searchEnd + 2000
        else:
            searchStart = searchStart - 2000
    if re.search("near-gene-3", snpInfo):
        if bedData.strand == '+':
            searchStart = searchStart - 500
        else:
            searchEnd = searchEnd + 500
    searchQuery = """select * from gafGeneXref where grch37LiteChrom = '%s'
                       and grch37LiteChromStart <= %d and grch37LiteChromEnd >= %d""" \
                       % (bedData.chrom, searchStart, searchEnd)
    cursor.execute(searchQuery)
    #
    # Be prepared that we might get multiple entries.  If we do, return
    # them separated with semicolons.
    #
    geneString = ""
    locusString = ""
    delimiter = ""
    for row in cursor.fetchall():
        if not re.search(re.sub("\?", "\?", row['geneName']), geneString):
            geneString = "%s%s%s" % (geneString, delimiter, row['geneName'])
            locusString = "%s%s%s" % (locusString, delimiter,
                                      row['grch37LiteLocus'])
            delimiter = ";"
    return((geneString, locusString))
    
    

db = MySQLdb.connect(host="localhost", db="hg19", read_default_file="~/.my.cnf")
cursor = db.cursor(MySQLdb.cursors.DictCursor)

parser = argparse.ArgumentParser()
parser.add_argument('grch37LiteBed', type=str, 
                    help="Input bed file in grch37-lite coordinates")
parser.add_argument('hg19Bed', type=str, 
                    help="Input bed file in hg19 coordinates")
parser.add_argument("-s", dest="snpTable", help="Table with the latest dbSNP data",
                    default="snp138")
parser.add_argument("-n", dest="entryNumber", help="Initial entry number",
                    default=0)
args = parser.parse_args()

entryNumber = args.entryNumber

#
# We assume that the hg19 coordinates and the GRCh37-lite coordinates
# are sorted in the same sort order: by chrom and then chromStart.
# There might be hg19 entries missing from the GRCh37 coordinates,
# as some pieces of hg19 are not in GRCh37-lite.  For each line in the
# GRCh37-lite file, read to the corresponding entry in the hg19 file.
# If we read to the end of the hg19 file without finding the
# entry we're looking for, print out a big error message.
hg19Fp = open(args.hg19Bed)
grch37LiteFp = open(args.grch37LiteBed)
for grch37LiteRow in grch37LiteFp:
    grch37LiteBed = Bed.Bed(grch37LiteRow.rstrip().split())
    for hg19Row in hg19Fp:
        hg19Bed = Bed.Bed(hg19Row.rstrip().split())
        if hg19Bed.name == grch37LiteBed.name:
            break
    if hg19Bed.name != grch37LiteBed.name:
        sys.exit("Error: missing entry for %s in the GRCh37-lite bed" % (hg19Bed.name))
    gg = Grch37LiteGaf.GafDbSnp(grch37LiteBed, entryNumber=entryNumber)
    gg.featureInfo = basicSnpInfo(hg19Bed, args.snpTable, cursor)
    (gg.gene, gg.geneLocus) = mapSnpToLocus(grch37LiteBed, gg.featureInfo, cursor)
    entryNumber = entryNumber + 1
    gg.entryNumber = entryNumber
    gg.write(sys.stdout)
exit(0)
