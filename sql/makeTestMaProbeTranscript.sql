CREATE TEMPORARY TABLE tmpUnchangedTranscripts AS
    SELECT kg.name FROM knownGene kg, knownGeneOld5 kgo
     WHERE kg.name = kgo.name AND kg.chrom = kgo.chrom
       AND kg.strand = kgo.strand AND kg.exonStarts = kgo.exonStarts
       AND kg.exonEnds = kgo.exonEnds
;
CREATE TEMPORARY TABLE gafMaProbeTranscriptData (
    probeName VARCHAR(100) NOT NULL,
    transcript VARCHAR(100) NOT NULL
);
LOAD DATA LOCAL INFILE "/cluster/home/cline/projects/TCGA-GAF-source/data/test/MAprobe.transcript.test.superset.txt" 
    INTO TABLE gafMaProbeTranscriptData 
    FIELDS TERMINATED BY ' '
;

SELECT g.* FROM gafMaProbeTranscriptData g, tmpUnchangedTranscripts t
  WHERE g.transcript = t.name
;
