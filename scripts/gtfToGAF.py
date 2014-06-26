#!/usr/bin/env python

import sys, os, re, getopt, argparse

from RangeFinder import RangeFinder
from GtfParse import FeatureObject, Gene, TranscriptTable, Transcript
import Grch37LiteGaf

parser = argparse.ArgumentParser(description="Create gene, transcript, exon level GAF files from input Gencode GTF file")
parser.add_argument('inputGtf', type=str,help="Input GTF file")
parser.add_argument('baseName', type=str,help="Base filename such as v4.0 or tmp. Outputs will be gene.genome.<baseName>.gaf, \
	transcript.genome.<baseName>.gaf and exon.genome.<baseName>.gaf")
parser.add_argument("-n", dest="entryNumber", help="Initial entry number",default="0")
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()

def gpGaf(gene, gFile, tFile, eFile, jFile, entryNumber):
    """Create GAF records for gene, transcripts, and exons from input gene object"""
    try:
	gene.getExons()
    except AttributeError:		# empty object: print nothing
	return entryNumber
    entryNumber += 1
    gg = Grch37LiteGaf.GafGene(gene, createFromGTF=True, entryNumber=entryNumber)
    gg.write(gFile)
    exonIds = set()
    junctionIds = set()
    for tx in gene.features:
        tx.getCDScoords()
        entryNumber += 1
        tg = Grch37LiteGaf.GafTranscript(tx, createFromGTF=True, entryNumber=entryNumber)
	tg.geneLocus = gg.geneLocus
	tg.gene = gg.gene
        tg.write(tFile)
	junction = False	# tracker for exon ends for creating junctions
        for e in tx.features:
       	    if e.descriptor == 'exon':
		if junction:
		    junction.finish(e.eId, e.start)
		    if junction.id not in junctionIds:
			# create GAF junction
                        entryNumber += 1
			jg = Grch37LiteGaf.GafJunction(junction, createFromJunction=True, entryNumber=entryNumber)
	                jg.geneLocus = gg.geneLocus
	                jg.gene = gg.gene
	                jg.write(jFile)
		    junctionIds.add(junction.id)
	        if e.eId not in exonIds:
                    entryNumber += 1
	            eg = Grch37LiteGaf.GafExon(e, createFromGTF=True, entryNumber=entryNumber)
	            eg.geneLocus = gg.geneLocus
	            eg.gene = gg.gene
	            eg.write(eFile)
	        exonIds.add(e.eId)
		junction = junctionInfo(e.eId, e.stop, e.chr, e.strand)
    return entryNumber

class junctionInfo(object):
    """Holder for junction information"""
    def __init__(self, id, pos, chr, strand):
	self.startId = id
	self.startPos = pos
	self.chr = chr
	self.strand = strand
    def finish(self, id, pos):
	self.endId = id
	self.endPos = pos
	self.id = (':').join([self.startId, self.endId])

# Main

gtfGenes = []	# list of gene IDs
gtfTranscripts = TranscriptTable()
curGene = Gene(Transcript=None, isEmpty = True)	# create empty gene object


gf=('.').join(['gene','genome',args.baseName,'gaf'])
tf=('.').join(['transcript','genome',args.baseName,'gaf'])
ef=('.').join(['exon','genome',args.baseName,'gaf'])
jf=('.').join(['junction','genome',args.baseName,'gaf'])
gFile = open(gf, 'w')
tFile = open(tf, 'w')
eFile = open(ef, 'w')
jFile = open(jf, 'w')
entryNumber = int(args.entryNumber)

f = open(args.inputGtf,'r')
for gtf_line in f:
    gtf_line = gtf_line.strip()
    if not gtf_line or gtf_line.startswith("#"):            # skip empty lines
        continue
    feats = FeatureObject(gtf_line)
    if feats.descriptor in ['transcript', 'gene']:	# skip non standard GTF lines
	continue
    if not 'basic' in feats.descDict['tag']:		# only keep basic gencode set
	continue
    if 'PAR' in feats.descDict['tag']:		# remove chrY pseudoautosomal region (not part of GRCh37-lite)
	continue
    if not gtfTranscripts.addFeatToTranscript(feats):
        # transcript does not yet exist in table: create and add first feature
        newTx = Transcript(feats)
	if not curGene.add(newTx):
	    entryNumber = gpGaf(curGene, gFile, tFile, eFile, jFile, entryNumber)
	    curGene = Gene(Transcript=newTx)
	    gtfTranscripts = TranscriptTable()	# reset table
	gtfTranscripts.add(newTx)
f.close()

# print last

entryNumber=gpGaf(curGene, gFile, tFile, eFile, jFile, entryNumber)
print "%d entries printed" % entryNumber
gFile.close()
tFile.close()
eFile.close()
jFile.close()
