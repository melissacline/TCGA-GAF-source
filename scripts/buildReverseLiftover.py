#!/usr/bin/env python2.7


import re
import sys

def printUsage():
    print "buildReverseLiftover.py - given a liftover chain, build the lift-back hain\n\
    usage:\n\
       buildReverseLiftover.py liftOver.chain > reverseLiftOver.chain\n\
    \n\
    where:\n\
       liftOver.chain maps the coordinates of genome A to genome B.\n\
       reverseLiftOver.chain maps the coordinates of genome B to genome A.\n"


if len(sys.argv) < 1:
    printUsage()
else:
    inputLiftOver = sys.argv[1]
    fp = open(inputLiftOver)
    for line in fp:
        line = line.rstrip()
        tokens = line.split()
        if len(tokens) > 3:
            # this is the chain declaration line.  We need to swap the
            # positions of the target and query genome data.
            (label, score, tName, tSize, tStrand, tStart, tEnd,
             qName, qSize, qStrand, qStart, qEnd, id) = tokens
            print label, score, qName, qSize, qStrand, qStart, qEnd, \
                  tName, tSize, tStrand, tStart, tEnd, id
        elif len(tokens) == 3:
            # This declares the size of a block and the lengths of any
            # gaps in target and query until the start of the next block.
            # We need to swap the positions of the target and query starts.
            (size, dt, dq) = tokens
            print size, dq, dt
        elif len(tokens) == 1:
            # This is the ungapped size of the last block.  All we need to
            # do is echo it back out.
            print tokens[0]
    
