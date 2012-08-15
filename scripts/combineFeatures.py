#!/usr/bin/env python

import argparse
import Gaf
import re
import sys



prevGaf = Gaf.Gaf()
prevFeatureId = ""
for line in sys.stdin.readlines():
    nextGaf = Gaf.Gaf(line)
    if nextGaf.featureId != prevFeatureId:
        if len(prevGaf.featureId) > 0:
            prevGaf.write(sys.stdout)
        prevGaf = nextGaf
        prevFeatureId = nextGaf.featureId
    else:
        prevGaf.combine(nextGaf)
if nextGaf.featureId == prevGaf.featureId:
    prevGaf.write(sys.stdout)

