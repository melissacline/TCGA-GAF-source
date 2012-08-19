CREATE TEMPORARY TABLE tmpSelectedSnps (name VARCHAR(100) NOT NULL);
LOAD DATA LOCAL INFILE "/cluster/home/cline/projects/TCGA-GAF-source/data/test/testDbSnp.txt" INTO TABLE tmpSelectedSnps;
select s.chrom, s.chromStart, s.chromEnd, s.name, s.score, s.strand, s.chromStart, s.chromEnd, "0", "1", s.chromEnd - s.chromStart, "0" 
  FROM snp135 s, tmpSelectedSnps t WHERE t.name = s.name
;
