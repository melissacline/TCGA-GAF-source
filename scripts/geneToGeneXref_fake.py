#!/usr/bin/python

import sys, os, re, getopt
from BasicStringStuff import *

usage = sys.argv[0]+""" <gene.genome gaf file>

Create geneName, grch37LiteChrom, grch37LiteChromStart and grch37LiteChromEnd
from input GAF, output a table ready for upload to gafGeneXref.
This program does NOT fill in hg19 coordinates, because those are not queried by the other programs.

"""


# Main
# read in command line and options
try:
    opts, args = getopt.getopt(sys.argv[1:], "dch")
except getopt.GetoptError:
    
        # print help information and exit:
    print usage
    print "ERROR did not recognize input\n"
    sys.exit(2)

for o, a  in opts:
    if o == "-h":
        print usage
        sys.exit()


if len(args) != 1:
    sys.exit(usage)

# Run program

f = open(args[0],'r')

for line in f:
    line = line.strip()
    fields = line.split("\t")
    geneName = fields[1]
    gafLocus = fields[16]
    chr,locus, strand = gafLocus.split(':')
    start, end = locus.split('-')
    start = int(start) -1
    t = "0"
    print ("\t").join([geneName,gafLocus,t,t,t,t,t,chr, str(start), end, strand,t ])
f.close()


