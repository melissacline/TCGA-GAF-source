DROP TABLE IF EXISTS gafGeneXref;
CREATE TABLE gafGeneXref(
    geneName		VARCHAR(255)	NOT NULL,
    grch37LiteLocus	VARCHAR(10000)	NOT NULL,
    clusterId		INTEGER,
    chrom		VARCHAR(200)	NOT NULL,
    chromStart  	INTEGER		NOT NULL,
    chromEnd    	INTEGER		NOT NULL,
    strand      	ENUM('+','-')	NOT NULL,
    alias               VARCHAR(200),
    KEY(clusterId),
    KEY(geneName),
    KEY(alias),
    KEY(chrom, chromStart, strand),
    KEY(chrom, chromEnd, strand)
);
