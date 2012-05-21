#!/usr/bin/env python


import re
import sys

def printUsage():
    print "buildLiftoverHg19ToGRCh37-lite.py - build a skeletal liftover for these genomes.\n\
    usage:\n\
       buildLiftoverHg19ToGRCh37-lite.py Grch37LiteReadme hg19ChromSizes\n\
    \n\
    where:\n\
       Grch37LiteReadme is the readme for GRCh37-lite and indicates which random contigs to include.\n\
       hg19ChromSizes is a two-column file with one line per chromosome or contig\n\
       in hg19, with the first column containing its name and the second \n\
       containing its size.\n\
   watchout:\n\
       This script builds a skeletal liftover file that assumes that each \n\
       chromosome or contig included is mapped without gaps.  There are some\n\
       gaps, and you will need to edit the output file to detail the gaps.\n"

def getRandoms(readMeFilename):
    """Return a hash of which random contigs to include, and the name to use in the lifted file"""
    fp = open(readMeFilename)
    targetRandoms = dict()
    for line in fp:
        line = line.rstrip()
        tokens = line.split()
        if len(tokens) == 2:
            if re.search("RANDOM_CTG", tokens[0]):
                value = tokens[1]
                key = re.split("\.", tokens[1])[0]
                key = key.lower()
                targetRandoms[key] = value
    fp.close()
    return(targetRandoms)

def printChainLine(referenceName, size, queryName, id):
    """Print a chain entry that assumes a full ungapped alignment for this chr"""
    print "chain 100", \
          referenceName, size, "+ 0", size, \
          queryName, size, "+ 0", size, id
    print size


if len(sys.argv) <= 1:
    printUsage()
else:
    Grch37LiteReadme = sys.argv[1]
    hg19ChromSizes = sys.argv[2]
    randomsToInclude = getRandoms(Grch37LiteReadme)
    #
    # Print out data for all the 'regular' chromosomes, the ones with
    # no underscores in the first token.  For the rest, if the first token
    # is in the target randoms list, then print out the line using the
    # name stored as a value in the target randoms dictionary.
    id = 1
    fp = open(hg19ChromSizes)
    for line in fp:
        line = line.rstrip()
        tokens = line.split()
        if not re.search("_", tokens[0]):
            printChainLine(tokens[0], tokens[1], tokens[0], id)
            id = id + 1
        else:
            subtokens = re.split("_", tokens[0])
            if randomsToInclude.has_key(subtokens[1]):
                printChainLine(tokens[0], tokens[1], 
                               randomsToInclude[subtokens[1]], id)
                id = id + 1

    
