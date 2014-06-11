#!/usr/bin/python

import sys, os, re, getopt

libFolder="/cluster/home/jeltje/lib/"
if libFolder not in sys.path: sys.path.append(libFolder)
from RangeFinder import RangeFinder
from GtfParse import GeneTable, FeatureObject, Gene, TranscriptTable, Transcript
import Grch37LiteGaf



usage = sys.argv[0]+""" <gtf format inputfile>

Create gene, transcript, exon level GAF files from input Gencode file

Options:         

"""

def gpGaf(gene, gFile, tFile, eFile):
    """Create GAF records for gene, transcripts, and exons from input gene object"""
    try:
	gene.getExons()
    except AttributeError:		# empty object: print nothing
	return False
    gg = Grch37LiteGaf.GafGene(gene, createFromGTF=True)
    gg.write(gFile)
    for tx in gene.features:
        tg = Grch37LiteGaf.GafTranscript(tx, createFromGTF=True)
	tg.geneLocus = gg.geneLocus
	tg.gene = gg.gene
        tg.write(tFile)
    for e in tx.features:
    	if e.descriptor == 'exon':
	    eg = Grch37LiteGaf.GafExon(e, createFromGTF=True)
	    eg.geneLocus = gg.geneLocus
	    eg.gene = gg.gene
	    eg.write(eFile)



# Main

# read in command line and options
try:
    opts, args = getopt.getopt(sys.argv[1:], "dc")
except getopt.GetoptError:
    
        # print help information and exit:
    print usage
    print "ERROR did not recognize input\n"
    sys.exit(2)


for o, a  in opts:
    if o == "-h":
        print usage
        sys.exit()

# Read in gtf line, create FeatureObject and if neccessary, a gene object, then add feature to gene

#print sys.path

gtfGenes = []	# list of gene IDs
gtfTranscripts = TranscriptTable()
curGene = Gene(Transcript=None, isEmpty = True)	# create empty gene object

if len(args) != 1:
    sys.exit(usage)
    
gFile = open('gene.genome.v4_0.gaf', 'w')
tFile = open('transcript.genome.v4_0.gaf', 'w')
eFile = open('exon.genome.v4_0.gaf', 'w')

f = open(args[0],'r')
for gtf_line in f:
    gtf_line = gtf_line.strip()
    if not gtf_line or gtf_line.startswith("#"):            # skip empty lines
        continue
    feats = FeatureObject(gtf_line)
    if feats.descriptor in ['transcript', 'gene']:	# skip non standard GTF lines
	continue
    if not 'basic' in feats.descDict['tag']:		# only keep basic gencode set
	continue
    if not gtfTranscripts.addFeatToTranscript(feats):
        # transcript does not yet exist in table: create and add first feature
        newTx = Transcript(feats)
	if not curGene.add(newTx):
	    gpGaf(curGene, gFile, tFile, eFile)
	    curGene = Gene(Transcript=newTx)
	    gtfTranscripts = TranscriptTable()	# reset table
	gtfTranscripts.add(newTx)
f.close()

# print last

gpGaf(curGene, gFile, tFile, eFile)
gFile.close()
tFile.close()
eFile.close()
