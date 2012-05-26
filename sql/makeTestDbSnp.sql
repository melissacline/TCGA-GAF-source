--
-- makeTestDbSnp.sql: select some test SNPs for comparing the GAF 2.1 
-- and GAF 3.0 data.  The GAF 2.1 data was built with dbSNP v131, while
-- the GAF 3.0 data was built with dbSNP v135.  Start with SNPs that exist
-- in both versions of dbSNP and have the same observed, class, and func
-- fields and coordinates.  Then, omit any with multiple coordinates.
--

--
-- For efficiency, start with a candidate pool of 1000 SNPs with the same
-- name and bin, both of which are indexed fields.  We'll cut that down later.
DROP TABLE IF EXISTS tmpTestSnps;
CREATE TEMPORARY TABLE tmpTestSnps AS
    SELECT snp135.name FROM snp131, snp135
     WHERE snp135.name = snp131.name
       AND snp135.bin = snp131.bin
      ORDER BY rand() limit 1000
;

--
-- Now select snps from this list for which the remaining features are the
-- same in both tables.
--
DROP TABLE IF EXISTS tmpBetterSnps;
CREATE TEMPORARY TABLE tmpBetterSnps AS
    SELECT t.name FROM tmpTestSnps t, snp131, snp135
     WHERE snp135.name = t.name AND snp131.name = t.name
       AND snp131.chrom = snp135.chrom
       AND snp131.chromStart = snp135.chromStart
       AND snp131.chromEnd = snp135.chromEnd
       AND snp131.observed = snp135.observed
       AND snp131.class = snp135.class
       AND snp131.func = snp135.func
;


--
-- Out of these SNPs, see how many genomic locations each of them has
-- (it should be sufficient to just look at snp135 for this), and select
-- the ones with just one location
--
DROP TABLE IF EXISTS tmpSnpLocationCounts;
CREATE TEMPORARY TABLE tmpSnpLocationCounts AS
    SELECT t.name, COUNT(*) AS loci
      FROM tmpBetterSnps t, snp135
     WHERE t.name = snp135.name
     GROUP BY t.name
;

SELECT name FROM tmpSnpLocationCounts
    WHERE loci = 1
;
