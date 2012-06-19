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
parser.add_argument('gaf', type=str, 
                    help="pre-miRNA to genome GAF file")
args = parser.parse_args()

db = MySQLdb.connect(host="localhost", db="hg19", user="hgcat",
                     passwd="S3attl3-S7u")
cursor = db.cursor(MySQLdb.cursors.DictCursor)

#
# First, make a dictionary of the hg19 coordinates of each pre-miRNA,
# indexed by the accession (name).
hg19Coordinates = dict()
hg19Fp = open(args.hg19Bed)
for line in hg19Fp:
    line = line.rstrip()
    bb = Bed.Bed(line.split())
    hg19Coordinates[bb.name] = bb
hg19Fp.close

#
# Next, read through each GAF entry.  At each entry, first look up
# the corresponding hg19 coordinates, which may be used later
# for comparison with other tables in hg19.  For each GAF entry,
# store the gene symbol and locus from the GAF file in the geneName
# and grch37LiteLocus columns.  Leave the cluster ID blank.  Store the
# feature ID in the alias column.
gafFp = open(args.gaf)
for line in gafFp:
    line = line.rstrip()
    preMiRnaGaf = Gaf.Gaf(line)
    (grch37LiteChrom, grch37LiteCoords,
     grch37LiteStrand) = preMiRnaGaf.compositeCoordinates.split(":")
    grch37LiteChromStart = grch37LiteCoords.split("-")[0]
    grch37LiteChromEnd = grch37LiteCoords.split("-")[-1]
    accession = preMiRnaGaf.featureId.split("|")[1]
    assert hg19Coordinates.has_key(accession)
    hg19Coords = hg19Coordinates[accession]
    tokens = bb.name.split(";")
    clusterId = tokens.pop()
    geneName = ";".join(tokens)
    locus = "%s:%d-%d:%s" % (bb.chrom, bb.chromStart + 1,
                             bb.chromEnd, bb.strand)
    cursor.execute("""INSERT INTO gafGeneXref (geneName, grch37LiteLocus, hg19Chrom,
                                               hg19ChromStart, hg19ChromEnd, hg19Strand,
                                               grch37LiteChrom, grch37LiteChromStart,
                                               grch37LiteChromEnd, grch37LiteStrand,
                                               alias)
                      VALUES ('%s', '%s',
                              '%s', %s, %s, '%s',
                              '%s', %s, %s, '%s',
                              '%s')""" \
                   % (preMiRnaGaf.gene, preMiRnaGaf.geneLocus, hg19Coords.chrom,
                      hg19Coords.chromStart, hg19Coords.chromEnd, hg19Coords.strand,
                      grch37LiteChrom, grch37LiteChromStart,
                      grch37LiteChromEnd, grch37LiteStrand,
                      preMiRnaGaf.featureId))
