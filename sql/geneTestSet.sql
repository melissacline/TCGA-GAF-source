create temporary table tmpExactTranscripts as
    select kgNew.name from knownGene kgNew, knownGeneOld5 kgOld
       where kgNew.name = kgOld.name
         and kgNew.txStart = kgOld.txStart and kgNew.txEnd = kgOld.txEnd
         and kgNew.cdsStart = kgOld.cdsStart and kgNew.cdsEnd = kgOld.cdsEnd
         and kgNew.exonStarts = kgOld.exonStarts 
         and kgNew.exonEnds = kgOld.exonEnds
;

create temporary table kgSummary as 
   select ki.* from knownIsoformsClustersMerged ki, tmpExactTranscripts et
    where ki.transcript = et.name
;

create temporary table tmpClusterTotals as 
    select clusterId, count(*) as total from knownIsoformsClustersMerged
     group by clusterId;
create index clusterId on tmpClusterTotals (clusterId);

create temporary table tmpExactCounts as 
    select clusterId, count(*) as exact from kgSummary
      group by clusterId;
create index clusterId on tmpExactCounts (clusterId);

create temporary table tmpSelectedClusters as 
    select t1.clusterId from tmpClusterTotals t1, tmpExactCounts t2
        where t1.clusterId = t2.clusterId and t1.total = t2.exact
;
create index clusterId on tmpSelectedClusters (clusterId);

create temporary table tmpSymbolsWithOneCluster as 
   select geneSymbol from 
      (select kx.geneSymbol, count(*) as clusters 
         from kgXref kx, knownIsoformsClustersMerged ki 
        where kx.kgID = ki.transcript 
         group by geneSymbol) as a  
   where clusters = 1;
create index geneSymbol on tmpSymbolsWithOneCluster (geneSymbol);

create temporary table gaf21Genes (geneName VARCHAR(100) PRIMARY KEY);
LOAD DATA LOCAL INFILE '/hive/users/cline/TCGA/GAF3.0/data/test/genes.in.gaf21.txt'
    INTO TABLE gaf21Genes;

create temporary table tmpTestGenes as 
    select distinct h.symbol, h.entrezId, t.clusterId
        from kgXref k, hgnc h, kgSummary kg, tmpSelectedClusters t, gafGeneXref g,
             gaf21Genes g21, tmpSymbolsWithOneCluster s1
        where h.symbol = k.geneSymbol and k.kgId = kg.transcript
              and kg.clusterId = t.clusterId 
              and g.clusterId = t.clusterId and g.geneName not like '%of%'
              and h.symbol = g21.geneName and h.symbol = s1.geneSymbol
        order by rand() limit 500
;

select symbol, entrezId 
    into outfile '/hive/users/cline/TCGA/GAF3.0/data/test/testGenes.txt'
    from tmpTestGenes order by symbol
;

create temporary table tmpTestTranscripts as
    select ki.transcript from tmpTestGenes tg, knownIsoformsClustersMerged ki
      where ki.clusterId = tg.clusterId
    order by rand() limit 500
;

select transcript
    into outfile '/hive/users/cline/TCGA/GAF3.0/data/test/testTranscripts.txt'
    from tmpTestTranscripts
    order by transcript
;
