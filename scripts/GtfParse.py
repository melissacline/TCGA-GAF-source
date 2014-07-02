""" GTF parser """

import sys
from RangeFinder import RangeFinder
from collections import defaultdict 

def overlap(x, y, a, b):
        """Takes four numbers representing ranges (in numerical order), returns False if the first number is smaller than the third or the second larger than the fourth"""
        if y< a or x>b:
                return False
        return True

def between(x, a, b):
        "Takes three numbers, returns True if the first number lies between the other two"
        if x<a or x >b:
                return False
        return True


def splitGtfId(gtfIds):
    """ extracts transcript and gene name from last column in gtf"""
    splitids = gtfIds.split('"')
    return splitids[1], splitids[3]

def getDescriptions(gtfIds):
    """ extracts key-value pairs from the last GTF field such as transcript and gene name"""
    # this appends items rather than overwrite them 
    gtfIds = gtfIds.rstrip(';')
    descDict = defaultdict(list)
    for i in gtfIds.split("; "):
    	s = filter(None,(i.split(" ")))		# removes empty entries
	v = s[1].replace('"','')		# remove quotes
	descDict[s[0]].append(v)
    return descDict

class FeatureObject(object):
    """ parses gtf line"""
    def __init__(self, gtf_line):
        feat = gtf_line.split('\t')
        self.chr = feat[0]
        self.defn =  feat[1]
        self.descriptor = feat [2]
        self.start= int(feat[3])
        self.stop= int(feat[4])
        self.score= feat[5]
        self.strand= feat[6]
        self.frame = feat[7]
	self.descDict = getDescriptions(feat[8])
	self.gId = self.descDict["gene_id"][0]
	self.tId = self.descDict["transcript_id"][0]
	self.geneType = self.descDict["gene_type"][0]
	self.geneStatus = self.descDict["gene_status"][0]
	self.geneSymbol = self.descDict["gene_name"][0]
	self.transcriptType = self.descDict["transcript_type"][0]
	self.transcriptStatus = self.descDict["transcript_status"][0]
	self.transcriptSymbol = self.descDict["transcript_name"][0]
	self.havanaGene = ''
	if "havana_gene" in self.descDict:
            self.havanaGene = self.descDict["havana_gene"][0]
	self.havanaTranscript = ''
	if "havana_transcript" in self.descDict:
            self.havanaTranscript = self.descDict["havana_transcript"][0]
	self.eId = ''
	if "exon_id" in self.descDict:
	    self.eId = self.descDict["exon_id"][0]
        self.inline = gtf_line

class Gene(object):
    def __init__(self, Transcript, isEmpty=False):
        """ initializes chromosome, strand, and Ids from transcript so we don't have to call these every time"""
	if isEmpty:
	    self.gId = "emptyGene"
	    return None
        self.features=[]
        self.exons = RangeFinder()
        self.cdsExons = RangeFinder()
        self.chr = Transcript.chr
        self.strand = Transcript.strand
	self.descDict = Transcript.descDict
	self.gId = self.descDict["gene_id"][0]
	self.tId = self.descDict["transcript_id"][0]
	self.geneType = self.descDict["gene_type"][0]
	self.geneStatus = self.descDict["gene_status"][0]
	self.geneSymbol = self.descDict["gene_name"][0]
	self.transcriptType = self.descDict["transcript_type"][0]
	self.transcriptStatus = self.descDict["transcript_status"][0]
	self.transcriptSymbol = self.descDict["transcript_name"][0]
	self.havanaGene = Transcript.havanaGene
	self.havanaTranscript = Transcript.havanaTranscript
        self.gId = Transcript.gId
        self.features.append(Transcript)
    def printout(self):
        for f in self.features:
                print f.inline
    def add(self, Transcript):
        """adds Transcript to Gene object"""
	if Transcript.gId == self.gId:
            self.features.append(Transcript)
	    return True
    def getExons(self):
	"""Get the subset of features that represent exons"""
	self.exonStarts = []
	self.exonEnds = []
	txStarts = []
	txEnds = []
        for t in self.features:
	    t.getExons()
	    self.exonStarts, self.exonEnds = flatten(self.exonStarts, self.exonEnds, t.exonStarts, t.exonEnds)
            assert(len(self.exonStarts) == len(self.exonEnds))
	    txStarts.append(t.locusStart)
	    txEnds.append(t.locusEnd)
	self.locusStart = min(txStarts)
	self.locusEnd = max(txEnds)

def flatten(gStarts, gEnds, tStarts, tEnds):
	"""Given two sets of exon ranges, extend the boundaries of the first set with those of the second and merge or add exons as appropriate"""
	removedTx = set()
	for i in xrange(len(tStarts)):
	    for k in xrange(len(gStarts)):
		if overlap(gStarts[k], gEnds[k], tStarts[i], tEnds[i]):
		    if i in removedTx:		# we're merging this exon and the previous one
			gStarts[k] = gStarts[k-1]
                        # change ends of all preceding exons that have the same start
                        c = 1
                        while c <= k:
                            if gStarts[k] == gStarts[k-c]:
                                gEnds[k-c] = max([gEnds[k],tEnds[i]])
                            c += 1
		    else:
		        gStarts[k]= min([gStarts[k],tStarts[i]])
		    gEnds[k]= max([gEnds[k],tEnds[i]])
		    removedTx.add(i)
	# add any tx exons that did not overlap any existing exons
	gStarts.extend(delete_by_indices(tStarts, removedTx))
	gEnds.extend(delete_by_indices(tEnds, removedTx))
	# Remove duplicate exons (happens when one tx exon overlaps multiple gene exons)
	return sorted(list(set(gStarts))), sorted(list(set(gEnds)))

def delete_by_indices(lst, indices):
    return [ lst[i] for i in xrange(len(lst)) if i not in set(indices) ]

class Transcript(object):
    def __init__(self, FeatureObject):
        """ initializes chromosome, strand, and Ids from FeatureObject so we don't have to call these every time"""
        self.features=[]
        self.exons = RangeFinder()
        self.cdsExons = RangeFinder()
        self.chr = FeatureObject.chr
        self.strand = FeatureObject.strand
	self.descDict = FeatureObject.descDict
	self.gId = self.descDict["gene_id"][0]
	self.tId = self.descDict["transcript_id"][0]
	self.geneType = self.descDict["gene_type"][0]
	self.geneStatus = self.descDict["gene_status"][0]
	self.geneSymbol = self.descDict["gene_name"][0]
	self.transcriptType = self.descDict["transcript_type"][0]
	self.transcriptStatus = self.descDict["transcript_status"][0]
	self.transcriptSymbol = self.descDict["transcript_name"][0]
	self.havanaGene = FeatureObject.havanaGene
	self.havanaTranscript = FeatureObject.havanaTranscript
        self.tId = FeatureObject.tId
        self.gId = FeatureObject.gId
        self.features.append(FeatureObject)
	self.hasStartCodon = True	# unless defined otherwise
	self.hasStopCodon = True
	self.firstFrame = None
    def addRefseq(self, refseqDict):
        self.refseqIds = refseqDict.get(self.tId) 	# this returns None if not possible
    def printout(self):
        for f in self.features:
                print f.inline
    def add(self, FeatureObject):
        """adds FeatureObject to Gene object"""
        self.features.append(FeatureObject)
    def getExons(self):
	"""Get the subset of features that represent exons"""
	self.exonStarts = []
	self.exonEnds = []
	self.eIds = []
	self.features =  sorted(self.features, key=lambda f: f.start) 	# sort features by start position
        for f in self.features:
            if (f.descriptor == "exon"):
                self.exonStarts.append(f.start)
                self.exonEnds.append(f.stop)
		self.eIds.append(f.eId) 
	self.locusStart = min(self.exonStarts)
	self.locusEnd = max(self.exonEnds)
	self.exonStarts = sorted(self.exonStarts)
	self.exonEnds = sorted(self.exonEnds)
    def makeCDSExons(self):
        """ take features list and create a list of exonstarts and stops only for the CDS (for determining overlap)
	also checks if the transcript has start_codon and stop_codon features"""
        self.cdsStarts = []
        self.cdsEnds= []
	hasStartCodon = False
	hasStopCodon = False
	lastFrame = None
        for f in self.features:
            if (f.descriptor == "CDS"):
		if self.firstFrame is None:
		    self.firstFrame = int(f.frame)
		lastFrame = int(f.frame)
                self.cdsStarts.append(f.start)
                self.cdsEnds.append(f.stop)
            elif (f.descriptor == "start_codon"):
                hasStartCodon = True
            elif (f.descriptor == "stop_codon"):
                hasStopCodon = True
	if f.strand == '-':
	    self.firstFrame = lastFrame
        if self.cdsStarts:	# if there are CDS exons in this transcript
	    if not hasStartCodon:
		self.hasStartCodon = False
	    if not hasStopCodon:
		self.hasStopCodon = False
	self.cdsStarts.sort()
	self.cdsEnds.sort()
    def getCDScoords(self):
        """Map CDS exons to transcript exons to find start and stop coordinate in transcript"""
        self.getExons()
        self.makeCDSExons()
	self.localStart = None
	self.localEnd = None
        if not self.cdsStarts:
            return False
        txCounter = 0
        cStart = self.cdsStarts[0]
        cEnd = self.cdsEnds[-1]
	if self.strand == '+':
            for a, b in zip(self.exonStarts, self.exonEnds):
                if between(cStart, a, b):
                    self.localStart = txCounter + (cStart - a + 1)
                if between(cEnd, a, b):
                    self.localEnd = txCounter + (cEnd - a + 1) + 3  	# add 3 for stop codon
                    break
                txCounter += b - a + 1
	elif self.strand == '-':
            for a, b in zip(reversed(self.exonStarts), reversed(self.exonEnds)):
                if between(cEnd, a, b):
                    self.localStart = txCounter + (b - cEnd + 1)
                if between(cStart, a, b):
                    self.localEnd = txCounter + (b - cStart + 1) +3	# add 3 for stop codon
                    break
                txCounter += b - a + 1

    def makeExons(self):
        """ take features list and create a list of exonstarts and stops. This means merging utr and exon features (for determining indentity)"""
        self.exonStarts = []
        self.exonEnds= []
        keepStart =''
        keepStop = 0
        for f in self.features:
            if (keepStart and (f.start < keepStart)):
                print >> sys.stderr, "ERROR, file must be sorted by start position, please run validate_gtf.pl -f on the input file."
                sys.exit(2)
            if (f.start - 1 > keepStop):         # if it is a new exon
                if(keepStop > 0):
                    self.exonStarts.append(keepStart)
                    self.exonEnds.append(keepStop)
                keepStart = f.start
            keepStop = f.stop
        # now add the last
        self.exonStarts.append(keepStart)
        self.exonEnds.append(keepStop)
                
class TranscriptTable(object):
    """ holds all transcript objects for a single gene"""
    def __init__(self):
        self.genes = []
    def add (self, gene):
        """ adds transcript object to transcriptTable"""
        self.genes.append(gene)
    def addFeatToTranscript(self, FeatureObject):
        """ tests if transcript is already in transcriptTable object, and adds to it"""
        for g in self.genes:
            if(g.tId == FeatureObject.tId):     # by definition, all features in a transcript have the same tId
                g.add(FeatureObject)
                return True
	

class GeneTable(object):
    """ holds all gene objects """
    def __init__(self):
        self.genes = []
    def add (self, gene):
        """ adds gene object to geneTable"""
    	self.genes.append(gene)
    def addTranscriptToGene(self, Transcript):
        """ tests if gene is already in GeneTable object, if so, adds to it"""
        for g in self.genes:
            if(g.gId == Transcript.gId):     # by definition, all features in a gene have the same gId
                g.add(Transcript)
                return True
