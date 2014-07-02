
import Gaf
from pycbio.hgdata import Bed
import re


class Grch37LiteGaf(Gaf.Gaf):
    """This class contains GAF objects in which the composite is the GRCh37-lite
    genome"""
    def __init__(self, inputLine, entryNumber=0, createFromBedInput=False, createFromGTF=True, createFromJunction=False):
        """Derive a Grch37LiteGaf object.  If there is an input bed line,
        then derive it from the bed object"""
        if not (createFromBedInput or createFromGTF):
            super(Grch37LiteGaf, self).__init__(line=inputLine,
                                                entryNumber=entryNumber)
        elif createFromBedInput:
            super(Grch37LiteGaf, self).__init__(line=None,
                                                entryNumber=entryNumber)
            self.compositeId = "GRCh37-lite"
            self.compositeType = "genome"
            self.compositeDbSource = "NCBI"
            self.compositeDbVersion = "GRCh37-lite"
            self.alignmentType = "pairwise"
            self._coordinatesFromBed(inputLine)
        elif createFromGTF or createFromJunction:
            super(Grch37LiteGaf, self).__init__(line=None, entryNumber=entryNumber)
            self.compositeId = "GRCh37-lite"
            self.compositeType = "genome"
            self.compositeDbSource = "NCBI"
            self.compositeDbVersion = "GRCh37-lite"
            self.alignmentType = "pairwise"
	    if not createFromJunction:
                self._getGafCoordsFromGtf(inputLine)
                self.featureDbSource = "Gencode"
                self.featureDbVersion = "V19"
                self.featureDbDate = '20130731'

    def _getCompositeCoordsFromGtf(self, GTF):
        """Get the CompositeCoordinates from a list of exonstarts and exonends"""
        genoCoords = []
	if not hasattr(GTF,'exonStarts'):	# in exons, exonstarts and exonends do not exist
	    GTF.exonStarts = [GTF.start,]
	    GTF.exonEnds = [GTF.stop,]
        for a, b in zip(GTF.exonStarts, GTF.exonEnds):
            genoCoords.append(('-').join([str(a),str(b)]))
            endpos = startpos + b - a
            txCoords.append(('-').join([str(startpos), str(endpos)]))
            startpos = endpos + 1
        self.compositeCoords = (':').join([GTF.chr, (',').join(i for i in genoCoords), GTF.strand])

    def _getGafCoordsFromGtf(self, GTF):
        """Get the CompositeCoordinates and map corresponding FeatureCoordinates from a list of exonstarts and exonends"""
        genoCoords = []
        txCoords = []
        startpos = 1
	if not hasattr(GTF,'exonStarts'):	# in exons, exonstarts and exonends do not exist
	    GTF.exonStarts = [GTF.start,]
	    GTF.exonEnds = [GTF.stop,]
        for a, b in zip(GTF.exonStarts, GTF.exonEnds):
            genoCoords.append(('-').join([str(a),str(b)]))
            endpos = startpos + b - a
            txCoords.append(('-').join([str(startpos), str(endpos)]))
            startpos = endpos + 1
        self.compositeCoordinates = (':').join([GTF.chr, (',').join(i for i in genoCoords), GTF.strand])
        if GTF.strand == '-':               # feature coords must be recalculated in reverse
            startpos = 1
            txCoords = []
            for a, b in zip(reversed(GTF.exonStarts), reversed(GTF.exonEnds)):
                endpos = startpos + b - a
                txCoords.append(('-').join([str(startpos), str(endpos)]))
                startpos = endpos + 1
        self.featureCoordinates = (',').join(i for i in txCoords)

    def _coordinatesFromBed(self, inputBed):
        """Given an input bed object, fill in the coordinate field of this
        GAF-derived object (and the name)"""
        self.featureId = inputBed.name
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
    def __init__(self, inputLine, entryNumber=0, createFromBedInput=False, createFromGTF=True):
        super(GafGene, self).__init__(inputLine, entryNumber=entryNumber,
                                            createFromBedInput=createFromBedInput, createFromGTF=createFromGTF)
        if createFromBedInput:
            self.featureType = "gene"
            self.featureDbSource = "calculated"
            self.featureSeqFileName = "genomic"
            self.gene = inputLine.name
            self.geneLocus = "%s:%d-%d:%s" % (inputLine.chrom, 
                                              inputLine.chromStart + 1,
                                              inputLine.chromEnd, inputLine.strand)
            self.featureInfo = "Confidence=%d" % (inputLine.score)
	elif createFromGTF:
            self.featureType = "gene"
            self.gene = ('|').join([inputLine.geneSymbol, inputLine.gId])
            self.geneLocus = "%s:%d-%d:%s" % (inputLine.chr, inputLine.exonStarts[0], inputLine.exonEnds[-1], inputLine.strand)
            self.featureInfo = "Gene_type=%s;Gene_status=%s;Gene_symbol=%s" % (inputLine.geneType, 
			inputLine.geneStatus, inputLine.geneSymbol)
            self.featureId = ('|').join([inputLine.geneSymbol, inputLine.gId])	# define here?
	    self.featureAliases = inputLine.havanaGene


class GafTranscript(Grch37LiteGaf):
    """GAF representation of a transcript on the GRCh37-lite genome"""
    def __init__(self, inputLine, entryNumber=0, createFromBedInput=False, createFromGTF=True):
        super(GafTranscript, self).__init__(inputLine, entryNumber=entryNumber,
                                            createFromBedInput=createFromBedInput, createFromGTF=createFromGTF)
        if createFromBedInput:
            tokens = inputBed.name.split(";")
            self.featureId = tokens.pop(0)
            self.gene = ";".join(tokens)
            self.featureType = 'transcript'
            self.featureDbSource = 'GencodeV19'
            self.featureDbDate = '20130731'
            self.featureSeqFileName = 'GencodeV19.fa'
            self.featureInfo = "Confidence=%d" % (inputLine.score)
            if inputLine.thickStart < inputLine.thickEnd:
                blockIdx = 0
                cdsStart = 0
                while (inputLine.blocks[blockIdx].end < inputLine.thickStart
                       and blockIdx < len(inputLine.blocks)):
                    cdsStart = cdsStart + inputLine.blocks[blockIdx].size
                    blockIdx = blockIdx + 1
                cdsStart = cdsStart + inputLine.thickStart - inputLine.blocks[blockIdx].start
                blockIdx = 0
                cdsStop = 0
                while inputLine.blocks[blockIdx].end < inputLine.thickEnd and blockIdx < len(inputLine.blocks):
                    cdsStop = cdsStop + inputLine.blocks[blockIdx].size
                    blockIdx = blockIdx + 1
                cdsStop = cdsStop + inputLine.thickEnd - inputLine.blocks[blockIdx].start
                if inputLine.strand == '-':
                    totalTranscriptSize = 0
                    for bb in inputLine.blocks:
                        totalTranscriptSize= totalTranscriptSize + bb.size
                    cdsStartCopy = cdsStart
                    cdsStart = totalTranscriptSize - cdsStop
                    cdsStop = totalTranscriptSize - cdsStartCopy
                self.featureInfo = "%s;CDSstart=%d;CDSstop=%d" \
                                   % (self.featureInfo, cdsStart + 1, cdsStop)
	elif createFromGTF:
            self.featureType = 'transcript'
       	    self.featureSeqFileName = 'GencodeV19.fa'       # contains all transcript sequences
            self.featureId = inputLine.tId
            self.featureInfo = "Transcript_type=%s;Transcript_status=%s;Transcript_symbol=%s" % (inputLine.transcriptType,
                                     inputLine.transcriptStatus, inputLine.transcriptSymbol)
            if inputLine.refseqIds:
                self.featureInfo += ";RefSeqId="+(";RefSeqId=").join(inputLine.refseqIds)
            if not inputLine.hasStartCodon:
                self.featureInfo += ";Warning=noStartCodon;FrameOffset=%d" % inputLine.firstFrame 
            if not inputLine.hasStopCodon:
                self.featureInfo += ";Warning=noStopCodon"
	    if inputLine.localStart:
                self.featureInfo += ";CDSstart=%d;CDSstop=%d" \
                                   % (inputLine.localStart, inputLine.localEnd)
	    self.featureAliases = inputLine.havanaTranscript


class GafExon(Grch37LiteGaf):
    """GAF representation of an exon on the GRCh37-lite genome"""
    def __init__(self, inputLine, entryNumber=0, createFromBedInput=False, createFromGTF=True,
                 exonType=""):
        super(GafExon, self).__init__(inputLine, entryNumber=entryNumber,
                                      createFromBedInput=createFromBedInput, 
				      createFromGTF=createFromGTF)
        if createFromBedInput:
            self.featureId = "%s:%d-%d:%s" \
                             % (inputLine.chrom, inputLine.chromStart + 1,
                                inputLine.chromEnd, inputLine.strand)
            tokens = inputBed.name.split(";")
            del tokens[0]
            self.gene = ";".join(tokens)
            self.featureType = exonType
            self.featureDbSource = "calculated"
            self.featureSeqFileName = "genomic"
	elif createFromGTF:
            self.featureType = 'exon'
            self.featureId = inputLine.eId

class GafJunction(Grch37LiteGaf):
    """GAF representation of a splice junction"""
    def __init__(self, inputLine, entryNumber=0, createFromBedInput=False, createFromJunction=True,
                 junction=0):
        super(GafJunction, self).__init__(inputLine, entryNumber=entryNumber,
                                      createFromBedInput=createFromBedInput, 
				      createFromJunction=createFromJunction)
        if createFromBedInput:
            self.featureType = "junction"
            self.featureDbSource = "calculated"
            self.featureSeqFileName = "genomic"
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
	elif createFromJunction:
            self.featureType = 'junction'
            self.featureDbSource = "calculated"
            self.featureCoordinates = "1,2"
	    self.featureId = inputLine.id
	    self.compositeCoordinates = "%s:%d,%d:%s" % (inputLine.chr, inputLine.startPos, 
					inputLine.endPos, inputLine.strand)
        
                



class GafPreMiRna(Grch37LiteGaf):
    """GAF representation of a pre-miRNA on the GRCh37-lite genome"""
    def __init__(self, inputBed, entryNumber=0, createFromBedInput=True):
        super(GafPreMiRna, self).__init__(inputBed, entryNumber=entryNumber,
                                          createFromBedInput=createFromBedInput)
        if createFromBedInput:
            self.featureType = "pre-miRNA"
            self.featureDbSource = "miRBase"
            self.featureDbVersion = "release19"
            self.featureDbDate = "20120812"
            self.featureSeqFileName = "miRBase19.hsa.gff3"
            self.alignmentType = "miRNA"

class GafMiRna(GafPreMiRna):
    """GAF representation of a miRNA on the GRCh37-lite genome"""
    def __init__(self, inputBed, entryNumber=0, createFromBedInput=True):
        super(GafMiRna, self).__init__(inputBed, entryNumber=entryNumber,
                                       createFromBedInput=createFromBedInput)
        if createFromBedInput:
            self.entryNumber=entryNumber
            self.featureType = "miRNA"

    def combine(self, newGafMiRna):
        """Combine a new GAF entry with this one, typically to handle one
        miRNAs with multiple genomic coordinates"""
        nextFeatureInfo = "%s,%s" % (self.featureInfo,
                                     re.sub("pre-miRNA=", "",
                                            newGafMiRna.featureInfo))
        super(GafMiRna, self).combine(newGafMiRna)
        self.featureInfo = nextFeatureInfo
        


class GafMaProbe(Grch37LiteGaf):
    """GAF representation of a microarray probe on the GRCh37-lite genome"""
    def __init__(self, inputBed, entryNumber=0, createFromBedInput=True):
        super(GafMaProbe, self).__init__(inputBed, entryNumber=entryNumber,
                                         createFromBedInput=createFromBedInput)
        if createFromBedInput:
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
    def __init__(self, inputBed, entryNumber=0, createFromBedInput=True):
        super(GafDbSnp, self).__init__(inputBed, entryNumber=entryNumber,
                                       createFromBedInput=createFromBedInput)
        if createFromBedInput:
            self.featureType = "dbSNP"
            self.featureDbSource = "UCSC"
            self.featureDbVersion = "v135"
            self.featureCoordinates = self.singleBaseCoordFixup(self.featureCoordinates)
            self.compositeCoordinates = self.singleBaseCoordFixup(self.compositeCoordinates)

    def _parseFeatureInfo(self, gafDbSnpRecord):
        """Given a dbSNP GAF record, parse out allele, class, and function
        from its FeatureInfo field and return them"""
        (alleleInfo, classInfo, functionInfo) = gafDbSnpRecord.featureInfo.split(";")
        assert re.search("^Alleles=", alleleInfo)
        assert re.search("^dbSNPclass=", classInfo)
        assert re.search("^dbSNPfunction=", functionInfo)
        alleleData = re.split("^Alleles=", alleleInfo)[1]
        classData = re.split("^dbSNPclass=", classInfo)[1]
        functionData = re.split("^dbSNPfunction=", functionInfo)[1]
        return((alleleData, classData, functionData))

    def combine(self, newGafDbSnp):
        """Combine a new GAF entry with this one, typically to handle one
        SNP with multiple genomic coordinates"""
        (myAlleles, myClass, myFunction) = self._parseFeatureInfo(self)
        (newAlleles, newClass, newFunction) = self._parseFeatureInfo(newGafDbSnp)
        nextFeatureInfo = "Alleles=%s,%s;dbSNPclass=%s,%s;dbSNPfunction=%s,%s" \
                          % (myAlleles, newAlleles, myClass, newClass,
                             myFunction, newFunction)
        super(GafDbSnp, self).combine(newGafDbSnp)
        self.featureInfo = nextFeatureInfo
        


class GafAffySnp(Grch37LiteGaf):
    def __init__(self, inputBed, entryNumber=0, createFromBedInput=True):
        super(GafAffySnp, self).__init__(inputBed, entryNumber=entryNumber,
                                         createFromBedInput=createFromBedInput)
        if createFromBedInput:
            self.featureType = "AffySNP"
            self.featureDbSource = "calculated"
            self.featureDbVersion = "GenomeWideSNP_6"
            self.featureDbSource = "Affymetrix"
            self.featureCoordinates = self.singleBaseCoordFixup(self.featureCoordinates)
            self.compositeCoordinates = self.singleBaseCoordFixup(self.compositeCoordinates)
        
            

