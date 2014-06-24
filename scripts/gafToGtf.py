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
        orf=dict(item.split("=") for item in featureInfo.split(";")) 
        try:
            self.cdsStart = int(orf['CDSstart'])
        except KeyError:
            self.cdsStart = None
        try:
            self.cdsStop = int(orf['CDSstop'])
        except KeyError:
            self.cdsStop = None
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
        bases of the first CDS). Also creates split starts and stops where necessary"""
        self.exonTable = []    # holder for exon objects
        # plus strand
        if not (self.cdsStart and self.cdsStop):
            return False
	if self.strand == '-':
            # if we're working on the minus strand, create the GTF lines in reverse
            cdsStart = copy.copy(self.cdsStart) -1
            cdsStop = copy.copy(self.cdsStop) -1
            curLabel = '5UTR'    # this label changes as we go through the gene
            curFrame = 0
            for e in reversed(self.exons):
                counter = 0
                gstart, gend = (int(i) for i in e.split('-'))
                gCdsStart = None
                gCdsStop = None
                for i in xrange(gend-gstart+1):
                    if counter == cdsStart:
                        gCdsStart = gend - counter # Found start, keep track
			print "found cdsStart", cdsStart, gCdsStart
                        cdsStart = -1
                    elif counter == cdsStop:
                        gCdsStop = gend - counter # Found end, keep track
			print "found cdsStop", cdsStop, gCdsStop
                        break    # we find the stop last, so stop trying
                    counter+=1
                # depending on whether the start and/or stop codon were found, create exons and assign labels
                if(gCdsStart and gCdsStop):    # CDS in single exon
                    if not gstart == gCdsStart:
                        self.exonTable.append(exon(gCdsStart+1, gend, '5UTR', '.'))
                    self.exonTable.append(exon(gCdsStart-2, gCdsStart, 'start_codon', '0'))
                    self.exonTable.append(exon(gCdsStop+3, gCdsStart, 'CDS', '0'))
                    self.exonTable.append(exon(gCdsStop, gCdsStop+2, 'stop_codon', '0'))
                    if not gend == gCdsStop+3:
                        self.exonTable.append(exon(gstart, gCdsStop-1, '3UTR', '.'))
                    curLabel = '3UTR'
                elif gCdsStart:    
                    if not gstart == gCdsStart:
                        self.exonTable.append(exon(gCdsStart+1, gend, '5UTR', '.'))
                    self.exonTable.append(exon(gCdsStart-2, gCdsStart, 'start_codon', '0'))
                    self.exonTable.append(exon(gstart, gCdsStart, 'CDS', '0'))
                    curLabel = 'CDS'
                    curFrame = (3 - ((gCdsStart -gstart +1 - curFrame) % 3)) %3
                elif gCdsStop:
                    self.exonTable.append(exon(gCdsStop+3, gend, 'CDS', curFrame))
                    self.exonTable.append(exon(gCdsStop, gCdsStop+2, 'stop_codon', '0'))
                    curLabel = '3UTR'
                    if not gstart == gCdsStop-1:
                        self.exonTable.append(exon(gstart, gCdsStop-1, curLabel, '.'))
                else:    # if nothing particular was found inside the exon, keep the previous label
                    self.exonTable.append(exon(gstart, gend, curLabel, curFrame))
                    curFrame = (3 - ((gend - gstart +1 - curFrame) % 3)) %3
                cdsStart -= (gend-gstart+1)
                cdsStop -= (gend-gstart+1)

############################################################
        elif self.strand == '+':
            cdsStart = copy.copy(self.cdsStart) -1
            cdsStop = copy.copy(self.cdsStop) -1
            curLabel = '5UTR'    # this label changes as we go through the gene
            curFrame = 0
            # Frame is calculated as (3 - ((length-frame) mod 3)) mod 3. 
            for e in self.exons:
                counter = 0
                gstart, gend = (int(i) for i in e.split('-'))
                gCdsStart = None
                gCdsStop = None
                # Here we go through the exon to see if the cdsStart or cdsStop are in it
                if cdsStop == 0:    # split stop; special case; set and deal with below
                    gCdsStop = gstart
                for i in xrange(gend-gstart+1):
                    if counter == cdsStart:
                        gCdsStart = gstart + counter # Found start, keep track
                        cdsStart = -1
                    elif counter == cdsStop:
                        gCdsStop = gstart + counter # Found end, keep track
                        break    # on plus strand we find the stop last, so stop trying
                    counter+=1
                # depending on whether the start and/or stop codon were found, create exons and assign labels
                if(gCdsStart and gCdsStop):    # CDS in single exon
                    if not gstart == gCdsStart:
                        self.exonTable.append(exon(gstart, gCdsStart-1, '5UTR', '.'))
                    self.exonTable.append(exon(gCdsStart, gCdsStart+2, 'start_codon', '0'))
                    self.exonTable.append(exon(gCdsStart, gCdsStop-3, 'CDS', '0'))
                    self.exonTable.append(exon(gCdsStop-2, gCdsStop, 'stop_codon', '0'))
                    if not gend == gCdsStop:
                        self.exonTable.append(exon(gCdsStop, gend, '3UTR', '.'))
                    curLabel = '3UTR'
                elif gCdsStart:    
                    if not gstart == gCdsStart:
                        self.exonTable.append(exon(gstart, gCdsStart-1, '5UTR', '.'))
                    if (gCdsStart+2 > gend):    # allow for split starts
                        curLabel = 'start_2'
                        curFrame = 1
                        if gCdsStart + 1 > gend:    # only first base is in this exon
                            curLabel = 'start_1'
                            curFrame = 2
                        self.exonTable.append(exon(gCdsStart, gend, 'start_codon', '0'))
                    else:
                        self.exonTable.append(exon(gCdsStart, gCdsStart+2, 'start_codon', '0'))
                        curLabel = 'CDS'
                    self.exonTable.append(exon(gCdsStart, gend, 'CDS', '0'))
                    curFrame = 0
                    # Frame for next exon is calculated as (3 - ((length-frame) mod 3)) mod 3. 
                    curFrame = (3 - ((gend - gCdsStart +1 - curFrame) % 3)) %3
                elif gCdsStop:
                    # if there's a split stop, the last added CDS exon must be changed
                    if (gCdsStop-1 < gstart):
                        e = self.exonTable.pop()
                        self.exonTable.append(exon(e.start, e.end-2, curLabel, e.readingFrame))
                        curLabel = 'stop_codon'
                        self.exonTable.append(exon(e.end-1, e.end, curLabel, '0'))
                        self.exonTable.append(exon(gCdsStop, gCdsStop, curLabel, '1'))
                        curLabel = '3UTR'
                        if not gCdsStop+1 == gend:
                            self.exonTable.append(exon(gCdsStop+1, gend, curLabel, '.'))
                    elif (gCdsStop-2 < gstart):
                        e = self.exonTable.pop()
                        self.exonTable.append(exon(e.start, e.end-1, curLabel, e.readingFrame))
                        curLabel = 'stop_codon'
                        self.exonTable.append(exon(e.end, e.end, curLabel, '0'))
                        self.exonTable.append(exon(gCdsStop-1, gCdsStop, curLabel, '2'))
                        curLabel = '3UTR'
                        if not gCdsStop+1 == gend:
                            self.exonTable.append(exon(gCdsStop+1, gend, curLabel, '.'))
                    else:
                        self.exonTable.append(exon(gstart, gCdsStop-3, 'CDS', curFrame))
                        self.exonTable.append(exon(gCdsStop-2, gCdsStop, 'stop_codon', '0'))
                    curLabel = '3UTR'
                elif curLabel == 'start_1':
                    self.exonTable.append(exon(gstart, gstart+1, 'start_codon', '2'))
                    curLabel = 'CDS'
                    curFrame = 2
                    self.exonTable.append(exon(gstart, gend, curLabel, curFrame))
                    curFrame = (3 - ((gend - gstart +1 - curFrame) % 3)) %3
                elif curLabel == 'start_2':
                    self.exonTable.append(exon(gstart, gstart, 'start_codon', '1'))
                    curLabel = 'CDS'
                    curFrame = 1
                    self.exonTable.append(exon(gstart, gend, curLabel, curFrame))
                    curFrame = (3 - ((gend - gstart +1 - curFrame) % 3)) %3
                else:    # if nothing particular was found inside the exon, keep the previous label
                    self.exonTable.append(exon(gstart, gend, curLabel, curFrame))
                    curFrame = (3 - ((gend - gstart +1 - curFrame) % 3)) %3
                cdsStart -= (gend-gstart+1)
                cdsStop -= (gend-gstart+1)
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

