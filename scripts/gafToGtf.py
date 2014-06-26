#!/usr/bin/python

import sys, os, re, getopt
cgiFolder="/hive/users/jeltje/lib"
if cgiFolder not in sys.path:
        sys.path.append(cgiFolder)
from BasicStringStuff import *

usage = sys.argv[0]+""" <file>

Take transcript to genome GAF file, output GTF with exons, 5UTR, start codon, stop codon, 3UTR.
Fun :-)

"""

class transcript(object):
    """Holds transcript objects"""
    def __init__(self, GAFline):
        """Extract transcript info from GAFline"""
        # use GAF field ID
        GAFline = GAFline.strip()
        fields=GAFline.split('\t')
        self.EntryNumber = int(fields[0])
        self.FeatureID = fields[1]
        CompositeCoordinates = fields[14]
        # Looks like chr21:45834475-45834827,45843259-45843364,45845033-45845446:-
        # split into chromosome, exons, and strand
        self.chr, ranges, self.strand = CompositeCoordinates.split(':')
        self.exons = ranges.split(',')
        FeatureCoordinates = fields[13]    # 1-551,552-570,571-1839
        self.txExons = FeatureCoordinates.split(',')
        self.Gene = fields[15]
        GeneLocus = fields[16].split(':')    # chr18:10454625-10488698:+
        self.geneStart, self.geneStop = (int(i) for i in GeneLocus[1].split('-'))
        featureInfo = fields[18]    # this may contain CDSstart CDSstop in format Confidence=400;CDSstart=353;CDSstop=607
	self._getWarnings(featureInfo)
        orf=dict(item.split("=") for item in featureInfo.split(";"))
        try:
            self.cdsStart = int(orf['CDSstart'])
        except KeyError:
            self.cdsStart = None
        try:
            self.cdsStop = int(orf['CDSstop'])
        except KeyError:
            self.cdsStop = None

    def _getWarnings(self,featureInfo):
	"""Warnings may occur multiple times in the featureInfo field. Get them out one by one"""
	self.hasNoStart = False
	self.hasNoStop = False
	for item in featureInfo.split(";"):
	    key, value = item.split("=") 
	    if key == "Warning":
	        if value == "noStartCodon":
		    self.hasNoStart = True
	        elif value == "noStopCodon":
		    self.hasNoStop = True
	    elif key == "FrameOffset":
		self.frameOffset = int(value)
    def makeExonGTF(self):
        """Creates exons in GTF format"""
        self.exonGTF = []
        # Must deal with duplicate gene names
        lastfield = "gene_id \"{0}\"; transcript_id \"{1}\";".format(self.Gene, self.FeatureID)
        head = ("\t").join([self.chr, 'GAF3', 'exon'])
        tail = ("\t").join(['.', self.strand, '.'])
        for e in self.exons:    
            start, end = e.split('-')
            outline = ("\t").join([head, start, end, tail, lastfield])
            self.exonGTF.append(outline)
    def mapCDS(self):
        """Splits exons in UTRs, CDS, and stop codon, and adds a start codon (first three 
        bases of the first CDS exon). Also creates split starts and stops where necessary"""
        self.exonTable = []    # holder for exon objects
        if not (self.cdsStart and self.cdsStop):
            return False
            # if we're working on the minus strand, create the GTF lines in reverse
        cdsStart = copy.copy(self.cdsStart) -1
        cdsStop = copy.copy(self.cdsStop) -1
        curLabel = '5UTR'    # this label changes as we go through the gene
        curFrame = 0
        splitStart = None
        splitStop = None
	if self.strand == '-':
            for e in reversed(self.exons):
                counter = 0
                gstart, gend = (int(i) for i in e.split('-'))
                gCdsStart = None
                gCdsStop = None
                # Here we go through the exon to see if the cdsStart or cdsStop are in it
                for i in xrange(gend-gstart+1):
                    if counter == cdsStart:
                        gCdsStart = gend - counter # Found start, keep track
			if gCdsStart - 2 < gstart:
                            splitStart = gCdsStart -gstart +1	# (number of start codon bases in current exon)
                        cdsStart = -1
                    elif counter == cdsStop:
                        gCdsStop = gend - counter # Found end, keep track
			if gCdsStop + 2 > gend:
                            splitStop = gend- gCdsStop +1 	# (number of stop codon bases in current exon)
                        break    # we find the stop last, so stop trying
                    counter+=1
		curFrame, curLabel= self.addGtfLinesMinStrand(gCdsStart, gCdsStop, 
                                        gstart, gend, curFrame, curLabel)
                cdsStart -= (gend-gstart+1)
                cdsStop -= (gend-gstart+1)
            # correct split starts and stop if found
            if splitStart:
                self.addSplitStartMinStrand(splitStart)
            if splitStop:
                self.addSplitStopMinStrand(splitStop)
        elif self.strand == '+':
            for e in self.exons:
                counter = 0
                gstart, gend = (int(i) for i in e.split('-'))
                gCdsStart = None
                gCdsStop = None
                for i in xrange(gend-gstart+1):
                    if counter == cdsStart:
                        gCdsStart = gstart + counter # Found start, keep track
			if gCdsStart + 2 > gend:
                            splitStart = gend -gCdsStart +1	# (number of start codon bases in current exon)
                        cdsStart = -1
                    elif counter == cdsStop:
                        gCdsStop = gstart + counter # Found end, keep track
			if gCdsStop - 2 < gstart:
                            splitStop = gCdsStop - gstart + 1 	# (number of stop codon bases in current exon)
                        break    # we find the stop last, so stop trying
                    counter+=1
		curFrame, curLabel= self.addGtfLinesPlusStrand(gCdsStart, gCdsStop, 
                                        gstart, gend, curFrame, curLabel)
                cdsStart -= (gend-gstart+1)
                cdsStop -= (gend-gstart+1)
            if splitStart:
                self.addSplitStartPlusStrand(splitStart)
            if splitStop:
                self.addSplitStopPlusStrand(splitStop)
        if (self.hasNoStart or self.hasNoStop):
            self.removeStartOrStop()

    def removeStartOrStop(self):
	"""Remove start codon or stop codon if there's a warning in the gaf that they were
	not present in the input GTF. Also check that there are no UTRs on that side"""
	removeMe = set()
	for i in xrange(len(self.exonTable)):
	    if (self.exonTable[i].label == "start_codon" and self.hasNoStart):
	        removeMe.add(i)
	    elif (self.exonTable[i].label == "stop_codon" and self.hasNoStop):
	        removeMe.add(i)
	    elif (self.exonTable[i].label == "5UTR" and self.hasNoStart):
	        print >>sys.stderr, "WARNING, found 5UTR in start codon free", self.FeatureID
	    elif (self.exonTable[i].label == "3UTR" and self.hasNoStop):
	        print >>sys.stderr, "WARNING, found 3UTR in stop codon free", self.FeatureID
	self.exonTable =  [ self.exonTable[i] for i in xrange(len(self.exonTable)) if i not in removeMe ]
	# correct frame in remaining CDS exons
	if self.hasNoStart:
	        curFrame = self.frameOffset
	        for i in self.exonTable:
		    if i.label == "CDS":
		        i.readingFrame = curFrame
                        curFrame = (3 - ((i.end - i.start +1 - curFrame) % 3)) %3


    def addGtfLinesPlusStrand(self, gCdsStart, gCdsStop, gstart, gend, curFrame, curLabel):
        """depending on whether the start and/or stop codon were found, create exons and assign labels"""
        if(gCdsStart and gCdsStop):    # CDS in single exon
            if gstart < gCdsStart-1:
                self.exonTable.append(exon(gstart, gCdsStart-1, '5UTR', '.'))
            self.exonTable.append(exon(gCdsStart, gCdsStart+2, 'start_codon', '0'))
            self.exonTable.append(exon(gCdsStart, gCdsStop-3, 'CDS', '0'))
            self.exonTable.append(exon(gCdsStop-2, gCdsStop, 'stop_codon', '0'))
            if not gend == gCdsStop:
                self.exonTable.append(exon(gCdsStop, gend, '3UTR', '.'))
            curLabel = '3UTR'
        elif gCdsStart:
            if gstart < gCdsStart-1:
                self.exonTable.append(exon(gstart, gCdsStart-1, '5UTR', '.'))
            self.exonTable.append(exon(gCdsStart, gCdsStart+2, 'start_codon', '0'))
            curLabel = 'CDS'
            self.exonTable.append(exon(gCdsStart, gend, 'CDS', '0'))
            curFrame = 0
            # Frame for next exon is calculated as (3 - ((length-frame) mod 3)) mod 3.
            curFrame = (3 - ((gend - gCdsStart +1 - curFrame) % 3)) %3
        elif gCdsStop:
            if (gstart <= gCdsStop-3):
                self.exonTable.append(exon(gstart, gCdsStop-3, 'CDS', curFrame))
            self.exonTable.append(exon(gCdsStop-2, gCdsStop, 'stop_codon', '0'))
	    if gCdsStop+1 < gend:
                self.exonTable.append(exon(gCdsStop+1, gend, '3UTR', '.'))
            curLabel = '3UTR'
        else:    # if nothing particular was found inside the exon, keep the previous label
            self.exonTable.append(exon(gstart, gend, curLabel, curFrame))
            curFrame = (3 - ((gend - gstart +1 - curFrame) % 3)) %3
        return curFrame, curLabel

    def addGtfLinesMinStrand(self, gCdsStart, gCdsStop, gstart, gend, curFrame, curLabel):
        """depending on whether the start and/or stop codon were found, create exons and assign labels"""
        if(gCdsStart and gCdsStop):    # CDS in single exon
            if gCdsStart + 1 < gend:
                self.exonTable.append(exon(gCdsStart+1, gend, '5UTR', '.'))
            self.exonTable.append(exon(gCdsStart-2, gCdsStart, 'start_codon', '0'))
            self.exonTable.append(exon(gCdsStop+3, gCdsStart, 'CDS', '0'))
            self.exonTable.append(exon(gCdsStop, gCdsStop+2, 'stop_codon', '0'))
            if gstart < gCdsStop-1:
                self.exonTable.append(exon(gstart, gCdsStop-1, '3UTR', '.'))
            curLabel = '3UTR'
        elif gCdsStart:
            if gCdsStart + 1 < gend:
                self.exonTable.append(exon(gCdsStart+1, gend, '5UTR', '.'))
            self.exonTable.append(exon(gCdsStart-2, gCdsStart, 'start_codon', '0'))
            self.exonTable.append(exon(gstart, gCdsStart, 'CDS', '0'))
            curLabel = 'CDS'
            curFrame = 0
            curFrame = (3 - ((gCdsStart -gstart +1 - curFrame) % 3)) %3
        elif gCdsStop:
            if gCdsStop + 3 <= gend:
                self.exonTable.append(exon(gCdsStop+3, gend, 'CDS', curFrame))
            self.exonTable.append(exon(gCdsStop, gCdsStop+2, 'stop_codon', '0'))
            curLabel = '3UTR'
            if gstart < gCdsStop-1:
                self.exonTable.append(exon(gstart, gCdsStop-1, curLabel, '.'))
        else:    # if nothing particular was found inside the exon, keep the previous label
            self.exonTable.append(exon(gstart, gend, curLabel, curFrame))
            curFrame = (3 - ((gend - gstart +1 - curFrame) % 3)) %3
        return curFrame, curLabel
    def addSplitStartPlusStrand(self, splitStart):
        """Correct exon table to add split start"""
	for i in xrange(len(self.exonTable)):
            if self.exonTable[i].label == 'start_codon':
		addStart = exon(self.exonTable[i+2].start,
                              self.exonTable[i+2].start + splitStart % 2, 'start_codon', 3-splitStart)
		self.exonTable.insert(i+2, addStart)
		return True
    def addSplitStopPlusStrand(self, splitStop):
        """Correct exon table to add split stop. The splitStop variable (1 or 2) holds information on split location"""
#	for i in xrange(len(self.exonTable)):
#	    print >>sys.stderr, self.exonTable[i].label, str(self.exonTable[i].start), str(self.exonTable[i].end)
	for i in xrange(len(self.exonTable)):
            if self.exonTable[i].label == 'stop_codon':
		self.exonTable[i].start = self.exonTable[i].end - (splitStop -1)
		self.exonTable[i].readingFrame = splitStop
		self.exonTable[i-1].end -= 3 - splitStop
		addStop = exon(self.exonTable[i-1].end + 1, self.exonTable[i-1].end + 1 + (splitStop %2 ), 
                      'stop_codon', 0)
#		del(self.exonTable[i-1])
		self.exonTable.insert(i, addStop)
		return True
    def addSplitStartMinStrand(self, splitStart):
        """Correct exon table to add split start"""
	for i in xrange(len(self.exonTable)):
            if self.exonTable[i].label == 'start_codon':
		self.exonTable[i].start = self.exonTable[i+1].start
		self.exonTable[i].end = self.exonTable[i+1].end
		addStart = exon(self.exonTable[i+2].end - splitStart % 2,
                              self.exonTable[i+2].end, 'start_codon', 3-splitStart)
		self.exonTable.insert(i+2, addStart)
		return True
    def addSplitStopMinStrand(self, splitStop):
        """Correct exon table to add split stop. The splitStop variable (1 or 2) holds information on split location"""
	for i in xrange(len(self.exonTable)):
            if self.exonTable[i].label == 'stop_codon':
                self.exonTable[i].readingFrame = splitStop
		self.exonTable[i-1].label = 'stop_codon'
                self.exonTable[i-1].readingFrame = 0
		self.exonTable[i-1].start = self.exonTable[i-2].start
		self.exonTable[i-1].end = self.exonTable[i-1].start + (splitStop % 2) # if one, add one, if two, add none
		self.exonTable[i-2].start = self.exonTable[i-1].end + 1 
		self.exonTable[i].end = self.exonTable[i].start + splitStop - 1
		return True

    def makeORFGTF(self):
        """Creates exons in GTF format"""
        self.orfGTF = []
        lastfield = "gene_id \"{0}\"; transcript_id \"{1}\";".format(self.Gene, self.FeatureID)
        head = ("\t").join([self.chr, 'GAF3'])
        tail = ("\t").join(['.', self.strand])
        for e in self.exonTable:    
            outline = ("\t").join([head, e.label, str(e.start), str(e.end), tail, str(e.readingFrame), lastfield])
            self.orfGTF.append(outline)

class exon(object):
    """Creates exon objects for GTF"""
    def __init__(self, start, end, label, readingFrame):
        self.start = start
        self.end = end
        self.label = label
        self.readingFrame = readingFrame
	if not self.label in ['CDS', 'start_codon', 'stop_codon']:
            self.readingFrame='.'


# Main
# read in command line and options
try:
    opts, args = getopt.getopt(sys.argv[1:], "dch")
except getopt.GetoptError:
    
        # print help information and exit:
    print usage
    print "ERROR did not recognize input\n"
    sys.exit(2)

for o, a  in opts:
#    if o == "-d":
#        doNotDelete = True
#    if o == "-c":
#        cdsOnly = True
    if o == "-h":
        print usage
        sys.exit()


if len(args) != 1:
    sys.exit(usage)

# Run program

f = open(args[0],'r')

for line in f:
    tx = transcript(line)
# must add check if fields are transcript and genome (or maybe otherwise just not produce CDS)
    tx.makeExonGTF()
    tx.mapCDS()
    tx.makeORFGTF()
    for l in tx.orfGTF:
        print l
    print
f.close()

