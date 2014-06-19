#!/usr/bin/env python

import argparse
from pycbio.hgdata import Bed
import Gaf
import MySQLdb
import MySQLdb.cursors
import sys

parser = argparse.ArgumentParser()
parser.add_argument('hg19Bed', type=str, 
                    help="Input bed file in hg19 coordinates")
parser.add_argument('grch37LiteBed', type=str, 
                    help="Input bed file in GRCh37-lite coordinates")
args = parser.parse_args()

db = MySQLdb.connect(host="localhost", db="hg19", read_default_file="~/.my.cnf")
cursor = db.cursor(MySQLdb.cursors.DictCursor)

#
# First, make a dictionary of the hg19 coordinates of each gene
hg19Coordinates = dict()
hg19Fp = open(args.hg19Bed)
for line in hg19Fp:
    line = line.rstrip()
    bb = Bed.Bed(line.split())
    hg19Coordinates[bb.name] = bb
hg19Fp.close

#
# Next, read through each entry in GRCh37-lite coordinates.
# At each entry, look up the hg19 coordinates, which will be
# used later for comparison with other tables in hg19.  Also at
# each entry, separate the gene name and cluster ID, build a
# GRCh37-lite locus string, and store the gene name, cluster ID,
# gene locus string, and hg19 coordinates within separate columns.
lineCount = 0
grch37LiteFp = open(args.grch37LiteBed)
for line in grch37LiteFp:
    line = line.rstrip()
    grch37LiteCoords = Bed.Bed(line.split())
    print "working on", grch37LiteCoords.name
    hg19Coords = hg19Coordinates[grch37LiteCoords.name]
    tokens = grch37LiteCoords.name.split(";")
    clusterId = tokens.pop()
    geneName = ";".join(tokens)
    locus = "%s:%d-%d:%s" % (grch37LiteCoords.chrom,
                             grch37LiteCoords.chromStart + 1,
                             grch37LiteCoords.chromEnd,
                             grch37LiteCoords.strand)
    cursor.execute("""INSERT INTO gafGeneXref (geneName, grch37LiteLocus, clusterId,
                                               hg19Chrom, hg19ChromStart,
                                               hg19ChromEnd, hg19Strand,
                                               grch37LiteChrom, grch37LiteChromStart,
                                               grch37LiteChromEnd, grch37LiteStrand)
                      VALUES ('%s', '%s', %s,
                              '%s', %s, %s, '%s',
                              '%s', %s, %s, '%s')""" \
                   % (geneName, locus, clusterId,
                      hg19Coords.chrom, str(int(hg19Coords.chromStart) + 1),
                      hg19Coords.chromEnd, hg19Coords.strand,
                      grch37LiteCoords.chrom, grch37LiteCoords.chromStart,
                      grch37LiteCoords.chromEnd, grch37LiteCoords.strand ))
    lineCount = lineCount + 1
    if lineCount % 100 == 0:
        print lineCount, "genes added"

#
# Merge data from entries with the same gene ID
cursor.execute("""SELECT geneName, count(*) AS loci FROM gafGeneXref
                  GROUP BY geneName""")
rowCount = 0
for row in cursor.fetchall():
    print "working on", row["geneName"], row["loci"], "entries"
    if row["loci"] > 1:
        basicGeneName = row["geneName"]
        totalCount = row["loci"]
        #
        # Build the target locus string
        locusString = ""
        delimiter = ""
        cursor.execute("SELECT grch37LiteLocus FROM gafGeneXref WHERE geneName = '%s'" \
                       % (basicGeneName))
        for locusRow in cursor.fetchall():
            locusString = locusString + delimiter + locusRow["grch37LiteLocus"]
            delimiter = ";"
        print "built locus string", locusString
        #
        # Rename the genes and update the locus string for all genes with
        # this gene name.  
        counter = 1
        cursor.execute("""SELECT clusterId FROM gafGeneXref 
                          WHERE geneName = '%s'""" % (basicGeneName))
        for idRow in cursor.fetchall():
            cursor.execute("""UPDATE gafGeneXref
                                 SET geneName = '%s|%dof%d', grch37LiteLocus = '%s'
                               WHERE clusterId = '%d'""" \
                           % (basicGeneName, counter, totalCount,
                              locusString, idRow["clusterId"]))
            counter = counter + 1
    rowCount = rowCount + 1
    if rowCount % 10 == 0:
        print rowCount, "loci verified/updated"
