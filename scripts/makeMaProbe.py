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
more than one.  From the old MAprobe GAF data, get the feature coordinate string.
While the bed files tell you what bases in the genome each probe is aligned to, the
old GAF data tells you which bases in the probe sequence aligned (in some cases, there
are gaps).  Generate a list with one maProbe-genome coordinate pair.  A subsequent
postprocessing step can combine composite coordinates for probes with more than one
alignment.

Inputs:
probes.grch37-lite.bed: alignments of the probes in GRCh37-lite coordinates
probes.hg19.bed: alignments of the probes in hg19 coordinates
oldMaProbe.gaf: the MAprobe-genome entries from the GAF 2.1 file

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
import Gaf
import MySQLdb
import MySQLdb.cursors
import sys

db = MySQLdb.connect(host="localhost", db="hg19", user="hgcat",
                     passwd="S3attl3-S7u")
cursor = db.cursor(MySQLdb.cursors.DictCursor)

parser = argparse.ArgumentParser()
parser.add_argument('grch37LiteBed', type=str,
                    help="Probe BED file, GRCh37-lite coordinates")
parser.add_argument('hg19Bed', type=str, help="Probe BED file, hg19 coordinates")
parser.add_argument('gaf21File', type=str, help="GAF2.1 file")
parser.add_argument("-n", dest="entryNumber", help="Initial entry number",
                    default=0)
args = parser.parse_args()

#
# Read the hg19 coordinates into a dictionary.  We'll use these
# for finding overlapping genes.
hg19Coords = dict()
hg19CoordsFp = open(args.hg19Bed, "rU")
for line in hg19CoordsFp:
    hg19Bed = Bed.Bed(line.rstrip().split())
    hg19Coords[hg19Bed.name] = hg19Bed
hg19CoordsFp.close()

#
# Next, read the old GAF file.  For each entry, store the feature coordinate
# string.  We'll use this information to tell us if there is a gap in the
# alignment of the feature.
featureCoords = dict()
gaf21Fp = open(args.gaf21File)
for line in gaf21Fp:
    gg21 = Gaf.Gaf()
    gg21.setFields(line.rstrip().split('\t'))
    featureCoords[gg21.featureId] = gg21.featureCoordinates
gaf21Fp.close()

#
# Finally, read and process each line in the GRCh37-lite coordinate file.
entryNumber = args.entryNumber
fp = open(args.grch37LiteBed)
for line in fp:
    line = line.rstrip()
    grch37LiteBed = Bed.Bed(line.split())
    hg19Bed = hg19Coords[grch37LiteBed.name]
    gg = Gaf.GafMaProbe(grch37LiteBed, entryNumber)

    #
    # Look up any overlapping genes.  If multiple overlapping genes are found,
    # parse the genes and gene locus strings into semicolon-delimited lists.
    geneXrefQuery = """SELECT geneName, grch37LiteLocus FROM gafGeneXref
                         WHERE chrom = '%s' AND chromStart <= %s
                           AND chromEnd >= %s
                           and strand = '%s'""" % (hg19Bed.chrom, hg19Bed.chromEnd,
                                                   hg19Bed.chromStart,
                                                   hg19Bed.strand)
    cursor.execute(geneXrefQuery)
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
