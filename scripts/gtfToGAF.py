#!/usr/bin/env python

import sys, os, re, getopt, argparse

from RangeFinder import RangeFinder
from GtfParse import FeatureObject, Gene, TranscriptTable, Transcript
import Grch37LiteGaf
import Gaf

parser = argparse.ArgumentParser(description="Create gene, transcript, exon level GAF files from input Gencode GTF file")
parser.add_argument('inputGtf', type=str,help="Input GTF file")
parser.add_argument('inputRefSeq', type=str,help="Input file with gencode and RefSeq IDs")
parser.add_argument('baseName', type=str,help="Base filename such as v4.0 or tmp. Outputs will be gene.genome.<baseName>.gaf, \
	transcript.genome.<baseName>.gaf and exon.genome.<baseName>.gaf")
parser.add_argument("-n", dest="entryNumber", help="Initial entry number",default="0")
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()

def gpGaf(gene, outputFiles, entryNumber, junctionNumber):
    """Create GAF records for gene, transcripts, and exons from input gene object"""
    try:
	gene.getExons()
    except AttributeError:		# empty object: print nothing
	return entryNumber, junctionNumber
    entryNumber += 1
    gg = Grch37LiteGaf.GafGene(gene, createFromGTF=True, entryNumber=entryNumber)
    gg.write(outputFiles.gFile)
    exonIds = set()
    junctionList = junctionInfo()
    for tx in gene.features:
        tx.getCDScoords()
        entryNumber += 1
        tg = Grch37LiteGaf.GafTranscript(tx, createFromGTF=True, entryNumber=entryNumber)
	tg.geneLocus = gg.geneLocus
	tg.gene = gg.gene
        tg.write(outputFiles.tFile)
        entryNumber += 1
        tgg = Gaf.Gaf(entryNumber=entryNumber)
        tgg.featureToComposite(tg, gg, fullBlocksOnly=False)
        tgg.write(outputFiles.tgFile)
	junct = False	# junction
        for e in tx.features:
       	    if e.descriptor == 'exon':
		if junct:
                    junct.finish(e.start)
                    junct = junctionList.add(junct, junctionNumber)	# if the junction is already there, this will return that junction
                    entryNumber += 1
                    j = Grch37LiteGaf.GafJunction(junct, createFromJunction=True, entryNumber=entryNumber)
                    jt = Gaf.Gaf(entryNumber=entryNumber)
                    jt.featureToComposite(j, tg, fullBlocksOnly=False)
                    jt.write(outputFiles.jtFile)
                entryNumber += 1
	        eg = Grch37LiteGaf.GafExon(e, createFromGTF=True, entryNumber=entryNumber)
	        eg.geneLocus = gg.geneLocus
	        eg.gene = gg.gene
	        if e.eId not in exonIds:
	            eg.write(outputFiles.eFile)
                    entryNumber += 1
                    egg = Gaf.Gaf(entryNumber=entryNumber)
                    egg.featureToComposite(eg, gg, fullBlocksOnly=False)
                    egg.write(outputFiles.egFile)
                    entryNumber += 1
                # output one gaf line per transcript
                et = Gaf.Gaf(entryNumber=entryNumber)
                et.featureToComposite(eg, tg, fullBlocksOnly=True)
                et.write(outputFiles.etFile)
	        exonIds.add(e.eId)
		junct = junction(e.chr, e.strand, e.stop)
    for j in junctionList.table:
        # create GAF junction
        entryNumber += 1
        jg = Grch37LiteGaf.GafJunction(j, createFromJunction=True, entryNumber=entryNumber)
        jg.geneLocus = gg.geneLocus
        jg.gene = gg.gene
        jg.write(outputFiles.jFile)
        entryNumber += 1
        jgg = Gaf.Gaf(entryNumber=entryNumber)
        jgg.featureToComposite(jg, gg, fullBlocksOnly=False)
        jgg.write(outputFiles.jgFile)
    return entryNumber, junctionNumber + len(junctionList.table)

class junction(object):
    """Holder for junction information"""
    def __init__(self, chr, strand, startPos):
	self.startPos = int(startPos)
	self.chr = chr
	self.strand = strand
    def finish(self, endPos):
	self.endPos = int(endPos)
    def name(self, idNr):
	self.id = "JUNC%06d" % idNr

class junctionInfo(object):
    """Holder for junction information"""
    def __init__(self):
	self.table = []
    def add(self, junct, number):
        """Test if junction is in table. If not, names the new junction and adds it """
        for s in self.table:
            if (junct.startPos == s.startPos and junct.endPos == s.endPos):
                return s
	junct.name(len(self.table) + number +1)
	self.table.append(junct)
	return junct

class OutFiles(object):
    def __init__(self, baseName):
        gf=('.').join(['gene','genome',baseName,'gaf'])
        tf=('.').join(['transcript','genome', baseName,'gaf'])
        tgf=('.').join(['transcript','gene', baseName,'gaf'])
        ef=('.').join(['exon','genome', baseName,'gaf'])
        etf=('.').join(['exon','transcript', baseName,'gaf'])
        egf=('.').join(['exon','gene', baseName,'gaf'])
        jf=('.').join(['junction','genome', baseName,'gaf'])
        jtf=('.').join(['junction','transcript', baseName,'gaf'])
        jgf=('.').join(['junction','gene', baseName,'gaf'])
        self.gFile = open(gf, 'w')
        self.tFile = open(tf, 'w')
        self.tgFile = open(tgf, 'w')
        self.eFile = open(ef, 'w')
        self.etFile = open(etf, 'w')
        self.egFile = open(egf, 'w')
        self.jFile = open(jf, 'w')
        self.jtFile = open(jtf, 'w')
        self.jgFile = open(jgf, 'w')
    def close(self):
        self.gFile.close()
        self.tFile.close()
        self.tgFile.close()
        self.eFile.close()
        self.etFile.close()
        self.egFile.close()
        self.jFile.close()
        self.jtFile.close()
        self.jgFile.close()

# Main

gtfGenes = []	# list of gene IDs
gtfTranscripts = TranscriptTable()
curGene = Gene(Transcript=None, isEmpty = True)	# create empty gene object

outputFiles = OutFiles(args.baseName)

entryNumber = int(args.entryNumber)
junctionNumber = 0
refseqDict = dict()
f = open(args.inputRefSeq,'r')
for rs_line in f:
    rs_line = rs_line.strip()
    fields = rs_line.split("\t")
    enst = fields[0]
    refseq = fields[1]
    if enst in refseqDict:
        refseqDict[enst].append(refseq)
    else:
        refseqDict[enst] = [refseq,]
f.close()

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
        newTx.addRefseq(refseqDict)
	if not curGene.add(newTx):
	    entryNumber, junctionNumber = gpGaf(curGene, outputFiles, entryNumber, junctionNumber)
	    curGene = Gene(Transcript=newTx)
	    gtfTranscripts = TranscriptTable()	# reset table
	gtfTranscripts.add(newTx)
f.close()

# print last

entryNumber, junctionNumber=gpGaf(curGene, outputFiles, entryNumber, junctionNumber)
print "%d entries printed" % entryNumber

outputFiles.close()
