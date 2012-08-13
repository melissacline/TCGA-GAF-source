#!/usr/bin/env python
"""
makeMiRna.py: generate the miRNA-genome GAF records

Usage:
makeMiRna.py -n 1000 miRna.genome.bed hsa.gff3 > miRNA.genome.gaf

Generate these records given the GRCh37-lite coordinates of the miRNA records
(which were obtained in hg19 coordinate spacefrom the miRBase hsa.gff3 file)
and the miRBase hsa.gff3 file (used for mapping between different ID types).

Each entry in the miRNA-genome BED file is named with the miRBase ID
(eg MIMAT0005890_1).  For each entry, use the key-value pairs to look up the
name (eg. hsa-miR-1302) and accession (eg MIMAT0005890).  The name and accession
will be concatenated to form the GAF ID (eg hsa-miR-1302|MIMAT0005890).  Also from
hsa.gff3, obtain the ID of the pre-miRNA that the miRNA is derived from
(eg MI0006369_1), and use the key-value pairs from hsa.gff3 to look up the name
and accession of the pre-miRNA.  List the pre-miRNA in the featureInfo column
(e.g. pre-miRNA=hsa-mir-1302-8|MI0006369)

Arguments:
miRnaBed: a BED file giving the miRNA GRCh37-lite coordinates
hsa.gff3: the coordinate and name data distributed by miRBase

Options
-n: Initial entry number (GFF column 1). This number starts at 1 by default
    and is incremented for each output line.

"""

import argparse
import Bio.SeqIO
from pycbio.hgdata import Bed
import Grch37LiteGaf
from BCBio import GFF
import MySQLdb
import MySQLdb.cursors
import re
import sys

def preMiRnaGeneLookup(preMiRnaAlias, cursor):
    """Given a pre-miRNA name (its alias in the gafGeneXref table), look up
       the gene name (HugoSymbol|EntrezGeneId) and the locus coordinate string"""
    query ="""SELECT geneName, grch37LiteLocus FROM gafGeneXref
               WHERE alias = '%s'""" % (preMiRnaAlias)
    cursor.execute(query)
    if cursor.rowcount != 1:
        print "rowcount", cursor.rowcount, "problem query", query
    assert(cursor.rowcount == 1)
    row = cursor.fetchone()
    return(row['geneName'], row['grch37LiteLocus'])


parser = argparse.ArgumentParser()
parser.add_argument('miRnaBed', type=str, help="Input bed file")
parser.add_argument('inputGff', type=str,  help="Input GFF file")
parser.add_argument("-n", dest="entryNumber", help="Initial entry number",
                    default=0)
args = parser.parse_args()

db = MySQLdb.connect(host="localhost", db="hg19", user="hgcat",
                     passwd="S3attl3-S7u")
cursor = db.cursor(MySQLdb.cursors.DictCursor)

entryNumber = args.entryNumber

#
# Read over the gff file collecting data for ID mapping. 
# map the ID to what will be the displayed ID: the combination of the
# name and accession number.  For miRNAs, also record the ID of the
# pre-miRNA that it was derived from.
miRnaToPreMiRna = dict()
idToLabel = dict()
gffIter = GFF.parse(args.inputGff)
for chrom in gffIter:
    for hit in chrom.features:
        id = hit.id
        label = "%s|%s" % (hit.qualifiers["Name"][0],
                           hit.qualifiers["accession_number"][0])
        idToLabel[id] = label
        if hit.type == "miRNA":
            miRnaToPreMiRna[hit.id] = hit.qualifiers["derives_from"][0]
            

#
# Read the bed file containing the GRCh37-lite coordinates.
# While converting each line to GAF format, replace the ID
# with the miRNA name, look up the name of the pre-miRNA that
# the miRNA is derived from, and note that in the featureInfo field.
miRnaBedFp = open(args.miRnaBed)
for line in miRnaBedFp:
    bb = Bed.Bed(line.rstrip().split())
    gg = Grch37LiteGaf.GafMiRna(bb)
    gg.featureId = idToLabel[gg.featureId]
    preMiRnaId = miRnaToPreMiRna[bb.name]
    preMiRnaName = idToLabel[preMiRnaId]
    gg.featureInfo = "pre-miRNA=%s" % (preMiRnaName)
    (geneName, geneLocus) = preMiRnaGeneLookup(preMiRnaName, cursor)
    gg.gene = geneName
    gg.geneLocus = geneLocus
    entryNumber = entryNumber + 1
    gg.entryNumber = entryNumber
    gg.write(sys.stdout)
exit(0)


