#!/usr/bin/env python

import argparse
from pycbio.hgdata import Bed
import sys

parser = argparse.ArgumentParser()
parser.add_argument('coordinate', type=int, help="Coordinate")
args = parser.parse_args()


for line in sys.stdin.readlines():
    bb = Bed.Bed(line.rstrip().split())
    relativePosition = 0
    if bb.strand == '+':
        for block in bb.blocks:
            if args.coordinate >= block.start and args.coordinate <= block.end:
                relativePosition = relativePosition + args.coordinate - block.start
                break
            else:
                relativePosition = relativePosition + block.end - block.start
    else:
        blockIdx = -1
        while bb.blocks[blockIdx].start > args.coordinate:
            blockLength = bb.blocks[blockIdx].end - bb.blocks[blockIdx].start
            relativePosition = relativePosition + blockLength
            blockIdx = blockIdx - 1
        lengthIntoThisBlock = bb.blocks[blockIdx].end - args.coordinate
        relativePosition = relativePosition + lengthIntoThisBlock
    print relativePosition
    
    
    
