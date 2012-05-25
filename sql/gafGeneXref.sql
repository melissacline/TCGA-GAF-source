DROP TABLE IF EXISTS gafGeneXref;
CREATE TABLE gafGeneXref(
    geneName		VARCHAR(255)	NOT NULL,
    grch37LiteLocus	VARCHAR(1000)	NOT NULL,
    clusterId		INTEGER,
    hg19Chrom		VARCHAR(200)	NOT NULL,
    hg19ChromStart  	INTEGER		NOT NULL,
    hg19ChromEnd    	INTEGER		NOT NULL,
    hg19Strand      	ENUM('+','-')	NOT NULL,
    grch37LiteChrom		VARCHAR(200)	NOT NULL,
    grch37LiteChromStart  	INTEGER		NOT NULL,
    grch37LiteChromEnd    	INTEGER		NOT NULL,
    grch37LiteStrand      	ENUM('+','-')	NOT NULL,
    alias               VARCHAR(200),
    KEY(clusterId),
    KEY(geneName),
    KEY(alias),
    KEY(hg19Chrom, hg19ChromStart, hg19Strand),
    KEY(hg19Chrom, hg19ChromEnd, hg19Strand),
    KEY(grch37LiteChrom, grch37LiteChromStart, grch37LiteStrand),
    KEY(grch37LiteChrom, grch37LiteChromEnd, grch37LiteStrand)
);
