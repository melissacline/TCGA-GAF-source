#!/usr/bin/env python
"""
gafToGff3: Given GAF input from stdin, translate it to GFF3 and send to stdout

Usage: inputfile.gaf | gafToGff3 > inputfile.gff3

Note: this script uses the following modules:
- the Gaf module (distributed along with it)
- the BioPython Seq, SeqFeature, and SeqRecord modules
- the BCBio GFF module (http://github.com/chapmanb/bcbb/tree/master/gff)

Assumptions:
"""

import Gaf
from BCBio import GFF
import re
import Bio.Seq
import Bio.SeqFeature
import Bio.SeqRecord
import sys



def buildQualifierList(gafRecord):
    """Given a GAF record, return a list of qualifiers for when
    the record will be turned into a SeqFeature object"""
    qualifiers = dict()
    if len(gafRecord.featureDbVersion) > 0:
        qualifiers["featureDbVersion"] = gafRecord.featureDbVersion
    if len(gafRecord.featureDbDate) > 0:
        qualifiers["featureDbDate"] = gafRecord.featureDbDate
    if len(gafRecord.featureSeqFileName) > 0:
        qualifiers["featureSeqFileName"] = gafRecord.featureSeqFileName
    if len(gafRecord.compositeDbSource) > 0:
        qualifiers["compositeDbSource"] = gafRecord.compositeDbSource
    if len(gafRecord.compositeDbVersion) > 0:
        qualifiers["compositeDbVersion"] = gafRecord.compositeDbVersion
    if len(gafRecord.compositeDbDate) > 0:
        qualifiers["compositeDbDate"] = gafRecord.compositeDbDate
    qualifiers["alignmentType"] = gafRecord.alignmentType
    qualifiers["gene"] = gafRecord.gene
    qualifiers["geneLocus"] = gafRecord.geneLocus
    if len(gafRecord.featureAliases) > 0:
        qualifiers["featureAliases"] = gafRecord.featureAliases
    featureInfoFields = gafRecord.featureInfo.split(";")
    for info in featureInfoFields:
        (key,value) = info.split("=")
        key = key.lower()
        qualifiers[key] = value
    return(qualifiers)



def segmentToSeqFeature(segment, gafRecord, strand):
    """Given a GAF record, a strand identifier, and a segment
    from the feature-coordinate alignment, return an appropriate
    SeqFeature object"""
    qualifiers = buildQualifierList(gafRecord)
    # Parse out the endpoints of the segment.  GAF coordiantes can
    # either be (start-end), or if start and end are the same, they
    # can be (start), with a single endpoint.
    endpoints = segment.split("-")
    if len(endpoints) == 1:
        start = end = int(endpoints[0])
    else:
        start = int(endpoints[0])
        end = int(endpoints[1])
    location=Bio.SeqFeature.FeatureLocation(Bio.SeqFeature.ExactPosition(start),
                                            Bio.SeqFeature.ExactPosition(end))
    sf = Bio.SeqFeature.SeqFeature(location=location, 
                                   type=gafRecord.featureType, strand=strand,
                                   id=gafRecord.featureId, qualifiers=qualifiers)
    return(sf)


def gafToSeqRecords(gafRecord):
    """Given this GAF record, return a SeqRecord iterator containing one
    or more SeqFeature objects.  There will be one SeqFeature object
    for each comma-delimited segment in the feature-composite alignment""" 
    if gafRecord.compositeType == "genome":
        # If this is a genomic alignment, then extract the chromosome and
        # strand from the composite coordinate string, and break the
        # composite coordinates into blocks.  Assert that the composite
        # coordinates string is properly-formed, with two colons delimiting
        # chrom:coordinates:strand.
        compositeCoordinateTokens = gafRecord.compositeCoordinates.split(":")
        assert(len(compositeCoordinateTokens) == 3)
        seqId = chromosome = compositeCoordinateTokens[0]
        if compositeCoordinateTokens[2] == "+":
            strand = 1
        else:
            assert(compositeCoordinateTokens[2] == "-")
            strand = -1
        compositeCoordinateSegments = compositeCoordinateTokens[1].split(",")
    else:
        # If this record is not a genomic alignment, make sure that the
        # composite string contains no colons (as it should not), 
        # set the SeqRecord seqId to the composite ID, and set the
        # strand to "." (undefined).
        assert(len(gafRecord.compositeCoordinates.split(":")) == 1)
        seqId = gafRecord.compositeId
        strand = 0
        compositeCoordinateSegments = gafRecord.compositeCoordinates.split(",")
    features = []
    for segment in compositeCoordinateSegments:
        sf = segmentToSeqFeature(segment, gafRecord, strand)
        features.append(sf)
    record = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(""), id=seqId,
                                     features=features)
    return(record)
                        

seqRecords = [] 
for line in sys.stdin.readlines():
    gafRecord = Gaf.Gaf(line)
    seqRecords.append(gafToSeqRecords(gafRecord))
GFF.write(seqRecords, sys.stdout)
