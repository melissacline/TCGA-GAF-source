DROP TABLE IF EXISTS gafGeneXref;
CREATE TABLE gafGeneXref(
    geneName	VARCHAR(255)	NOT NULL,
    locus	VARCHAR(10000)	NOT NULL,
    clusterId	INTEGER		NOT NULL,
    chrom	VARCHAR(200)	NOT NULL,
    chromStart  INTEGER		NOT NULL,
    chromEnd    INTEGER		NOT NULL,
    strand      ENUM('+','-')	NOT NULL,
    KEY(clusterId),
    KEY(geneName),
    KEY(chrom, chromStart, strand),
    KEY(chrom, chromEnd, strand)
);
