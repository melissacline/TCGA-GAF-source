
import Gaf
from pycbio.hgdata import Bed
import re


class Grch37LiteGaf(Gaf.Gaf):
    """This class contains GAF objects in which the composite is the GRCh37-lite
    genome"""
    def __init__(self, inputBed):
        """Derive a Grch37LiteGaf object from a Bed object"""
        super(Grch37LiteGaf, self).__init__()
        self.compositeId = "GRCh37-lite"
        self.compositeType = "genome"
        self.compositeDbSource = "NCBI"
        self.compositeDbVersion = "GRCh37-lite"
        self.alignmentType = "pairwise"
        delimiter = ""
        compositeCoord = ""
        if inputBed.blocks is not None:
            for block in inputBed.blocks:
                compositeCoord = compositeCoord + delimiter \
                                 + str(block.start + 1) \
                                 + "-" + str(block.end)
                delimiter = ","
        compositeCoord = inputBed.chrom + ":" + compositeCoord \
                         + ":" + inputBed.strand
        featureCoord = ""
        delimiter = ""
        thisBlockStart = 0
        if inputBed.blocks is not None:
            if inputBed.strand == '+':
                for block in inputBed.blocks:
                    featureCoord = featureCoord + delimiter \
                                   + str(thisBlockStart + 1) \
                                   + "-" + str(thisBlockStart + block.size)
                    delimiter = ","
                    thisBlockStart = thisBlockStart + block.size
            else:
                for ii in range(len(inputBed.blocks)-1, -1, -1):
                    block = inputBed.blocks[ii]
                    featureCoord = featureCoord + delimiter \
                                   + str(thisBlockStart + 1) \
                                   + "-" + str(thisBlockStart + block.size)
                    delimiter = ","
                    thisBlockStart = thisBlockStart + block.size
        self.featureCoordinates = featureCoord
        self.compositeCoordinates = compositeCoord
            
            
class GafGene(Grch37LiteGaf):
    """GAF representation of a gene on the GRCh37-lite genome"""
    def __init__(self, inputBed, entryNumber=0):
        super(GafGene, self).__init__(inputBed)
        self.entryNumber = entryNumber
        self.featureId = inputBed.name
        self.featureType = "gene"
        self.featureDbSource = "calculated"
        self.featureSeqFileName = "genomic"
        self.alignmentType = "pairwise"
        self.gene = inputBed.name
        self.geneLocus = inputBed.chrom + ":" + str(inputBed.chromStart+1) \
                         + "-" + str(inputBed.chromEnd) + ":" + inputBed.strand
        self.featureInfo = "Confidence=%d" % (inputBed.score)

class GafTranscript(Grch37LiteGaf):
    """GAF representation of a transcript on the GRCh37-lite genome"""
    def __init__(self, inputBed, entryNumber=0):
        super(GafTranscript, self).__init__(inputBed)
        self.entryNumber = entryNumber
        tokens = inputBed.name.split(";")
        self.featureId = tokens.pop(0)
        self.gene = ";".join(tokens)
        self.featureType = 'transcript'
        self.featureDbSource = 'UCSCgene'
        self.featureDbDate = '20120131'
        self.featureSeqFileName = 'UCSCgene.Jan2012.fa'
        self.featureInfo = "Confidence=%d" % (inputBed.score)
        if inputBed.thickStart < inputBed.thickEnd:
            blockIdx = 0
            cdsStart = 0
            while inputBed.blocks[blockIdx].end < inputBed.thickStart and blockIdx < len(inputBed.blocks):
                cdsStart = cdsStart + inputBed.blocks[blockIdx].size
                blockIdx = blockIdx + 1
            cdsStart = cdsStart + inputBed.thickStart - inputBed.blocks[blockIdx].start
            blockIdx = 0
            cdsStop = 0
            while inputBed.blocks[blockIdx].end < inputBed.thickEnd and blockIdx < len(inputBed.blocks):
                cdsStop = cdsStop + inputBed.blocks[blockIdx].size
                blockIdx = blockIdx + 1
            cdsStop = cdsStop + inputBed.thickEnd - inputBed.blocks[blockIdx].start
            if inputBed.strand == '-':
                totalTranscriptSize = 0
                for bb in inputBed.blocks:
                    totalTranscriptSize= totalTranscriptSize + bb.size
                cdsStartCopy = cdsStart
                cdsStart = totalTranscriptSize - cdsStop
                cdsStop = totalTranscriptSize - cdsStartCopy
            self.featureInfo = "%s;CDSstart=%d;CDSstop=%d" \
                               % (self.featureInfo, cdsStart + 1, cdsStop)

class GafExon(Grch37LiteGaf):
    """GAF representation of an exon on the GRCh37-lite genome"""
    def __init__(self, inputBed, entryNumber=0, exonType=""):
        super(GafExon, self).__init__(inputBed)
        self.entryNumber = entryNumber
        self.featureId = "%s:%d-%d:%s" \
                           % (inputBed.chrom, inputBed.chromStart + 1,
                              inputBed.chromEnd, inputBed.strand)
        tokens = inputBed.name.split(";")
        del tokens[0]
        self.gene = ";".join(tokens)
        self.featureType = exonType
        self.featureDbSource = "calculated"
        self.featureSeqFileName = "genomic"
        self.alignmentType = "pairwise"

class GafPreMiRna(Grch37LiteGaf):
    """GAF representation of a pre-miRNA on the GRCh37-lite genome"""
    def __init__(self, inputBed, entryNumber=0):
        super(GafPreMiRna, self).__init__(inputBed)
        self.entryNumber = entryNumber
        self.featureId = inputBed.name 
        self.featureType = "pre-miRNA"
        self.featureDbSource = "miRBase"
        self.featureDbVersion = "release18"
        self.featureDbDate = "20111108"
        self.featureSeqFileName = "miRBase18.hsa.gff3"
        self.alignmentType = "miRNA"

class GafMiRna(GafPreMiRna):
    """GAF representation of a miRNA on the GRCh37-lite genome"""
    def __init__(self, inputBed, entryNumber=0):
        super(GafMiRna, self).__init__(inputBed)
        self.entryNumber=entryNumber
        self.featureType = "miRNA"

class GafMaProbe(Grch37LiteGaf):
    """GAF representation of a microarray probe on the GRCh37-lite genome"""
    def __init__(self, inputBed, entryNumber=0):
        super(GafMaProbe, self).__init__(inputBed)
        self.entryNumber = entryNumber
        self.featureId = inputBed.name
        self.featureType = "MAprobe"
        self.featureDbSource = "AgilentG4502A_07_3"
        self.featureSeqFileName = "unc.edu_AgilentG4502A_07_3.fa"
        self.alignmentType = "pairwise"


class GafDbSnp(Grch37LiteGaf):
    """ Given a BED representation of a SNP, return its GAF
    representation.  There is one special case for SNPs.  If any block
    in the feature or composite coordiante string has size 1, with the
    start and end coordinate being the same, then only give that
    coordinate once.  In other words, instead of having feature
    coordinate 1-1 and composite coordinate chr15:1000-1000:+, set
    featureCoordinate to 1 and composite coordinate to chr15:1000:+
    """
    def __init__(self, inputBed, entryNumber=0):
        super(GafDbSnp, self).__init__(inputBed)
        self.entryNumber = entryNumber
        self.featureId = inputBed.name
        self.featureType = "dbSNP"
        self.featureDbSource = "UCSC"
        self.featureDbVersion = "v135"
        self.alignmentType = "pairwise"
        self.featureCoordinates = self.singleBaseCoordFixup(self.featureCoordinates)
        self.compositeCoordinates = self.singleBaseCoordFixup(self.compositeCoordinates)


class GafAffySnp(Grch37LiteGaf):
    def __init__(self, inputBed, entryNumber=0):
        super(GafAffySnp, self).__init__(inputBed)
        self.entryNumber = entryNumber
        self.featureId = inputBed.name
        self.featureType = "AffySNP"
        self.featureDbSource = "calculated"
        self.featureDbVersion = "GenomeWideSNP_6"
        self.alignmentType = "pairwise"
        self.featureDbSource = "Affymetrix"
        self.featureCoordinates = self.singleBaseCoordFixup(self.featureCoordinates)
        self.compositeCoordinates = self.singleBaseCoordFixup(self.compositeCoordinates)
        
            

class GafJunction(Grch37LiteGaf):
    """GAF representation of a splice junction"""
    def __init__(self, inputBed, entryNumber=0, junction=0):
        super(GafJunction, self).__init__(inputBed)
        self.entryNumber = entryNumber
        self.featureType = "junction"
        self.featureDbSource = "calculated"
        self.featureSeqFileName = "genomic"
        self.alignmentType = "pairwise"
        self.featureCoordinates = "1,2"
        self.featureId = "%s:%d:%s,%s:%d:%s" \
                         % (inputBed.chrom, inputBed.blocks[junction].end,
                            inputBed.strand, inputBed.chrom,
                            inputBed.blocks[junction+1].start+1,
                            inputBed.strand)
        self.compositeCoordinates = "%s:%d,%d:%s" \
                                    % (inputBed.chrom,
                                       inputBed.blocks[junction].end,
                                       inputBed.blocks[junction+1].start+1,
                                       inputBed.strand)

        
                

