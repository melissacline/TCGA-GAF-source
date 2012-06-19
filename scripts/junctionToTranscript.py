#!/usr/bin/env python

import argparse
import Gaf
import sys

    
def hasCommonEndpoints(junction, transcript):
    """Verify that the endpoints of the junction are endpoints of blocks
    in the transcript"""
    jBlocks = junction.gafToBlocks()
    assert len(jBlocks) == 2
    junctionStart = jBlocks[0].compositeEnd
    junctionEnd = jBlocks[1].compositeStart
    tBlocks = transcript.gafToBlocks()
    for ii in range(len(tBlocks) - 1):
        if tBlocks[ii].compositeEnd == junctionStart \
           and tBlocks[ii+1].compositeStart == junctionEnd:
            return(True)
    return(False)
    

parser = argparse.ArgumentParser()
parser.add_argument('junctionGaf', type=str, 
                    help="GAF containing junctions mapped against the genome")
parser.add_argument('transcriptGaf', type=str, 
                    help="GAF containing transcripts mapped against the genome")
parser.add_argument("-n", dest="entryNumber", help="Initial entry number",
                    default=0)
args = parser.parse_args()

entryNumber = args.entryNumber
#
# Read through the transcripts.  Make a dictionary with a list of transcripts
# for each gene.
geneToTranscript = dict()
transcriptFp = open(args.transcriptGaf)
for line in transcriptFp:
    transcript = Gaf.Gaf(line)
    if geneToTranscript.has_key(transcript.gene):
        geneToTranscript[transcript.gene].append(transcript)
    else:
        geneToTranscript[transcript.gene] = [transcript]
transcriptFp.close()

junctionFp = open(args.junctionGaf)
for line in junctionFp:
    junction = Gaf.Gaf(line)
    if geneToTranscript.has_key(junction.gene):
        transcriptsThisGene = geneToTranscript[junction.gene]
        for transcript in transcriptsThisGene:
            if hasCommonEndpoints(junction, transcript):
                junctionToTranscript = Gaf.Gaf()
                junctionToTranscript.featureToComposite(junction, transcript)
                assert len(junctionToTranscript.featureCoordinates) > 0
                junctionToTranscript.entryNumber = entryNumber
                print junctionToTranscript
                entryNumber = entryNumber + 1
