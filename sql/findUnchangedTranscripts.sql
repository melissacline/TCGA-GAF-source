select k5to6.newId from kg5ToKg6 k5to6, knownGene kg6, knownGeneOld5 kg5
     where k5to6.status = 'exact' 
       and k5to6.newId = kg6.name and k5to6.oldId = kg5.name 
       and kg5.cdsStart = kg6.cdsStart and kg5.cdsEnd = kg6.cdsEnd
       and kg5.txStart = kg6.txStart and kg5.txEnd = kg6.txEnd
;
