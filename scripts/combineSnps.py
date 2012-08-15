#!/usr/bin/env python

import argparse
import Grch37LiteGaf
import re
import sys


parser = argparse.ArgumentParser()
args = parser.parse_args()

prevSnpId = ""
for line in sys.stdin.readlines():
    nextGaf = Grch37LiteGaf.GafDbSnp(line.rstrip(), createFromBedInput=False)
    if nextGaf.featureId != prevSnpId:
        if len(prevGaf.featureId) > 0:
            prevGaf.write(sys.stdout)
        prevGaf = nextGaf
        prevSnpId = nextGaf.featureId
    else:
        prevGaf.combine(nextGaf)
if nextGaf.featureId != prevGaf.featureId:
    nextGaf.write(sys.stdout)
