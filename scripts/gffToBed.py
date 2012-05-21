#!/usr/bin/env python
"""
gffToBed.py: given a GFF file, output its contents in BED format.

Translate GFF to BED coordinates.  Optionally output only those entries
of a type or feature (GFF field 3) that matches the type parameter.
Use the type as the BED name by default, or optionally specify an element from
the attributes field (which are semicolon-delimited key-value pairs) to use as
the name.

Arguments:
inputGff: name of the GFF file

Options:
-t <type>: specifies a type of GFF record to translate.  The default is to
           translate everything.

-n <name>: specifies an attribute from the attributes field to use as the
           name.  The type is used as the name by default.  If a name
           attribute is specified but is not found in some record, there
           will be an assertion failure.

Usage:
gffToBed.py -t miRNA -n ID hsa.gff > miRNA.bed

Assumptions:

This code assumes that each line in the GFF file represents a single
ungapped alignment, equivalent to a single block in BED format.  This
does not support generating a multi-block BED from GFF.

Note that the type and name parameters are both case-sensitive.

"""

import argparse
from BCBio import GFF
import re

parser = argparse.ArgumentParser()
parser.add_argument('inputGff', type=str,  help="Input GFF file")
parser.add_argument('-t', dest="type", default=".*", help="Type to output")
parser.add_argument('-n', dest="name", default="",
                    help="Name of the field to use in the BED name field")
args = parser.parse_args()
typeSearchString = "^%s$" % (args.type)

gffIter = GFF.parse(args.inputGff)
for chrom in gffIter:
    chromName = "chr" + chrom.id
    for hit in chrom.features:
        if re.search(typeSearchString, hit.type):
            chromStart = hit.location.nofuzzy_start
            chromEnd = hit.location.nofuzzy_end
            if hit.strand:
                strand = "+"
            else:
                strand = "-"
            if args.name == "":
                hitName = hit.type
            else:
                hitName = hit.qualifiers[args.name][0]
            print "%s\t%d\t%d\t%s\t1\t%s\t%d\t%d\t0\t1\t%d\t0" \
                  % (chromName, int(chromStart), int(chromEnd),  hitName,
                     strand, int(chromStart), chromEnd,
                     int(chromEnd) - int(chromStart))
