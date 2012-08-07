#!/usr/bin/env python
"""
makeMaProbe.py: make the MAprobe-genome GAF records

Usage:
makeMaProbe.py -n 1000 probes.grch37-lite.bed probes.hg19.bed oldMaProbe.gaf \
    > maProbe.genome.gaf

Description:
Assemble the data for the MAprobe:genome records.  Get the coordinates from the
GRCH37-lite bed file.  From the hg19 coordinates, query the gafGeneXref table to
get any overlapping genes (assembled in a semicolon-delimited list, if there are
more than one.  Generate a list with one maProbe-genome coordinate pair.  A subsequent
postprocessing step can combine composite coordinates for probes with more than one
alignment.

Inputs:
probes.grch37-lite.bed: alignments of the probes in GRCh37-lite coordinates
probes.hg19.bed: alignments of the probes in hg19 coordinates

Options:
-n: starting record number (which by default starts at 1 and increments with each
    output record.

Hidden surprises:
In cases where a block of the composite alignment includes just a single base, this
is transformed to a start-end pair where start and end are the same coordinate.  For
example: chr2:1000,1220-1259:+ would become chr2:1000-1000,1220-1259:+.  Other types
of GAF records don't follow this convention.
"""

import argparse
from pycbio.hgdata import Bed
import Grch37LiteGaf
import MySQLdb
import MySQLdb.cursors
import sys

db = MySQLdb.connect(host="localhost", db="hg19", user="hgcat",
                     passwd="S3attl3-S7u")
cursor = db.cursor(MySQLdb.cursors.DictCursor)

parser = argparse.ArgumentParser()
parser.add_argument('grch37LiteBed', type=str,
                    help="Probe BED file, GRCh37-lite coordinates")
parser.add_argument("-n", dest="entryNumber", help="Initial entry number",
                    default=0)
parser.add_argument("-d", dest="debug", help="Optional debugging info", default=False)
args = parser.parse_args()



entryNumber = args.entryNumber
fp = open(args.grch37LiteBed)
for line in fp:
    line = line.rstrip()
    grch37LiteBed = Bed.Bed(line.split())
    gg = Grch37LiteGaf.GafMaProbe(grch37LiteBed, entryNumber)

    #
    # Look up any overlapping genes.  If multiple overlapping genes are found,
    # parse the genes and gene locus strings into semicolon-delimited lists.
    geneXrefQuery = """SELECT DISTINCT geneName, grch37LiteLocus
                          FROM gafGeneXref
                         WHERE grch37LiteChrom = '%s' AND grch37LiteChromStart <= %s
                           AND grch37LiteChromEnd >= %s
                           and grch37LiteStrand = '%s'""" % (grch37LiteBed.chrom,
                                                             grch37LiteBed.chromEnd,
                                                             grch37LiteBed.chromStart,
                                                             grch37LiteBed.strand)
    cursor.execute(geneXrefQuery)
    if (args.debug):
        print "executing query", geneXrefQuery
    delimiter  = ""
    gg.gene = ""
    gg.geneLocus = ""
    for row in cursor.fetchall():
        gg.gene = "%s%s%s" % (gg.gene, delimiter, row["geneName"])
        gg.geneLocus = "%s%s%s" % (gg.geneLocus, delimiter, row["grch37LiteLocus"])
        delimiter = ";"
    entryNumber = entryNumber + 1
    gg.entryNumber = entryNumber
    gg.write(sys.stdout)
exit(0)
