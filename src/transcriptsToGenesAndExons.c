/* transcriptsToGenesAndExons - given a set of UCSC Genes transcripts, generate
 * the component and composite exons, and the overall gene models. */
#include "common.h"
#include "bed.h"
#include "dystring.h"
#include "hash.h"
#include "hdb.h"
#include "hgConfig.h"
#include "linefile.h"
#include "options.h"
#include "rangeTree.h"

static struct optionSpec options[] = {
    {NULL, 0},
};

void usage()
{
    errAbort(
	"transcriptsToGenesAndExons - given a database table of transcripts,\n"
	"  generate the exon and gene model data\n"
        "usage:\n"
        "   transcriptsToGenesAndExons database transcriptTable clusterTable componentExons compositeExons transcripts genes geneXref\n"
        "where\n"
        "   database - the browser database with the appropriate gene table\n"
        "   transcriptTable - the transcript track to build from\n"
        "   clusterTable - a table which clusters transcripts by gene\n"
        "   componentExons - an output file for the component exon data\n"
        "   compositeExons - an output file for the composite exon data\n"
	"   transcripts - an output file for the transcript data\n"
        "   genes - an ouptut file for the gene models\n"
    );
}

struct hash *transcriptsThisChromAndStrand(struct sqlConnection *conn, char *transcripts,
    char *chrom, char strand)
    /* Get the set of transcripts for this chrom/strand */
{
    struct sqlResult *sr;
    char extraBuf[256], *extra, **row;
    int rowOffset;
    struct hash *geneStructHash = hashNew(16);
    struct genePred *gp;

    safef(extraBuf, sizeof(extraBuf), "strand = '%c'", strand);
    extra = extraBuf;
    sr = hChromQuery(conn, transcripts, chrom, extra, &rowOffset);
    while ((row = sqlNextRow(sr)) != NULL) {
	gp = genePredLoad(row + rowOffset);
	hashAdd(geneStructHash, gp->name, (void *) gp);
    }
    sqlFreeResult(&sr);
    return(geneStructHash);
}


char *clusterIdToGeneSymbol(struct sqlConnection *conn, int clusterId)
/* Given a cluster ID, return the gene symbol of the canonical transcript */
{
    struct dyString *geneSymbolQuery; 
    struct sqlResult *geneSymbolResult;
    char *geneSymbol, **row;

    geneSymbol = NULL;
    geneSymbolQuery = dyStringNew(256);
    dyStringPrintf(geneSymbolQuery, "select kx.geneSymbol from knownCanonical kc, kgXref kx where kc.transcript = kx.kgID and kc.clusterId = '%d'", clusterId);
    geneSymbolResult = sqlMustGetResult(conn, geneSymbolQuery->string);
    row = sqlNextRow(geneSymbolResult);
    geneSymbol = replaceChars(row[0], " ", "_");
    sqlFreeResult(&geneSymbolResult);
    dyStringFree(&geneSymbolQuery);
    return(geneSymbol);
}

int isSymbolInHugo(struct sqlConnection *conn, char *geneSymbol)
/* Return TRUE or FALSE depending on whether the symbol is in the hgnc table */
{
    struct dyString *hgncQuery;
    struct sqlResult *hgncResult;
    char **row;
    int returnValue;

    hgncQuery = dyStringNew(256);
    dyStringPrintf(hgncQuery, "select * from hgnc where symbol='%s'", 
		   geneSymbol);
    hgncResult = sqlGetResult(conn, hgncQuery->string);
    if ((row = sqlNextRow(hgncResult)) != NULL) {
	returnValue = TRUE;
    } else {
	returnValue = FALSE;
    }
    sqlFreeResult(&hgncResult);
    dyStringFree(&hgncQuery);
    return(returnValue);
}

char *geneSymbolFromSynonymOrPrvSymbol(struct sqlConnection *conn, char *geneSymbol)
/* Whatever is in geneSymbol is not a symbol in the hgnc table.  See
 * if it hits anything in the synonym or previous symbol column */
{
    struct dyString *symbolCountQuery, *symbolLookupQuery; 
    struct sqlResult *symbolLookupResult;
    char *newGeneSymbol, **row;

    newGeneSymbol = NULL;
    symbolCountQuery = dyStringNew(256);
    dyStringPrintf(symbolCountQuery, 
		   "hgnc where synonyms rlike '(%%,)*%s(,%%)*' or prvSymbols rlike '(%%,)*%s(,%%)*'", 
		   geneSymbol, geneSymbol);
    if (sqlRowCount(conn, symbolCountQuery->string) == 1) {
	symbolLookupQuery = dyStringNew(256);
	dyStringPrintf(symbolLookupQuery, 
		       "select symbol, entrezId from hgnc where synonyms rlike '(%%,)*%s(,%%)*' or prvSymbols rlike '(%%,)*%s(,%%)*'", 
		       geneSymbol, geneSymbol);
	symbolLookupResult = sqlGetResult(conn, symbolLookupQuery->string);
	row = sqlNextRow(symbolLookupResult);
	assert(row != NULL);
	newGeneSymbol = cloneString(row[0]);
	sqlFreeResult(&symbolLookupResult);
    }
    return(newGeneSymbol);
}



char *geneSymbolToEntrezId(struct sqlConnection *conn, char *geneSymbol)
/* Look up the entrez gene ID from the hgnc table, using the contents of 
 * geneSymbol to access the symbol field.  If this works, then we have the
 * entrez ID AND the HUGO-approved gene symbol */
{
    struct dyString *entrezIdQuery; 
    struct sqlResult *entrezIdResult;
    char *entrezGeneId, **row;

    entrezGeneId = NULL;
    entrezIdQuery = dyStringNew(256);
    dyStringPrintf(entrezIdQuery, "select entrezId from hgnc where symbol='%s'", 
		   geneSymbol);
    entrezIdResult = sqlGetResult(conn, entrezIdQuery->string);
    if ((row = sqlNextRow(entrezIdResult)) != NULL) {
	if (strcmp(row[0], "") != 0) {
	    entrezGeneId = cloneString(row[0]);
	}
    }
    sqlFreeResult(&entrezIdResult);
    dyStringFree(&entrezIdQuery);
    return(entrezGeneId);
}

char *clusterIdToEntrezId(struct sqlConnection *conn, int clusterId)
/* Given a cluster ID, look for the  LocusLink (Entrez Gene) ID for the cluster */
{
    struct dyString *locusLinkQuery; 
    struct sqlResult *locusLinkResult;
    char *entrezId, **row;

    entrezId = NULL;
    locusLinkQuery = dyStringNew(256);
    dyStringPrintf(locusLinkQuery, "select kl.value from knownCanonical kc, knownToLocusLink kl where kc.clusterId = %d and kc.transcript = kl.name", clusterId);
    locusLinkResult = sqlGetResult(conn, locusLinkQuery->string);
    if ((row = sqlNextRow(locusLinkResult)) != NULL) {
	entrezId = cloneString(row[0]);
    }
    sqlFreeResult(&locusLinkResult);
    dyStringFree(&locusLinkQuery);
    return(entrezId);
}


char *entrezIdToGeneSymbol(struct sqlConnection *conn, char *entrezId)
/* Look up the gene symbol from the hgnc table, using the entrez gene ID */
{
    struct dyString *geneSymbolQuery; 
    struct sqlResult *geneSymbolResult;
    char *newGeneSymbol, **row;

    newGeneSymbol = NULL;
    geneSymbolQuery = dyStringNew(256);
    dyStringPrintf(geneSymbolQuery, "select symbol from hgnc where entrezId='%s'", 
		   entrezId);
    geneSymbolResult = sqlGetResult(conn, geneSymbolQuery->string);
    if ((row = sqlNextRow(geneSymbolResult)) != NULL) {
	newGeneSymbol = cloneString(row[0]);
    }
    sqlFreeResult(&geneSymbolResult);
    dyStringFree(&geneSymbolQuery);
    return(newGeneSymbol);
}


char *clusterIdToGeneSymbolViaEnsembl(struct sqlConnection *conn, 
				    char *knownIsoformTable, 
				    int clusterId)
/* Look up the symbol from hgnc via knownToEnsembl */
{
    struct dyString *geneSymbolQuery;
    struct sqlResult *geneSymbolResult;
    char *geneSymbol, **row;

    /* Given the cluster ID, look up the transcripts from the knownIsoforms table.
     * Join the transcript IDs with those in the knownToEnsembl table, producing
     * a list of ENSEMBL transcripts.  Use those to get the ENSEMBL gene accession
     * from the ensGene table.  Use the gene accession to query into the hgnc table
     * to produce a symbol and an Entrez ID.  
     * Taking this route, there can be multiple mappings between cluster ID and 
     * symbol|entrezID.  Because of this, select only the one most frequent mapping.*/
    geneSymbol = NULL;
    geneSymbolQuery = dyStringNew(256);
    dyStringPrintf(geneSymbolQuery, "select h.symbol, count(*) as matches from %s ki, knownToEnsembl ke, ensGene eg, hgnc h where ki.transcript = ke.name and ke.value = eg.name and eg.name2 = h.ensId and ki.clusterId = %d group by h.symbol, h.entrezId order by matches limit 1", knownIsoformTable, clusterId);
    geneSymbolResult = sqlGetResult(conn, geneSymbolQuery->string);
    if ((row = sqlNextRow(geneSymbolResult)) != NULL) {
	if (strcmp(row[0], "") != 0) {
	    geneSymbol = cloneString(row[0]);
	} 
    }
    dyStringFree(&geneSymbolQuery);
    sqlFreeResult(&geneSymbolResult);
    return(geneSymbol);
}


char *geneNameThisCluster(struct sqlConnection *conn, int clusterId, 
			  char *knownIsoformTable)
{
    struct dyString *symbolBuf = dyStringNew(80);
    char *geneSymbol, *newGeneSymbol, *geneName;
    char *entrezGeneId = NULL;
    int hasHugoSymbol = FALSE;

    /* Start with the simple case.  Given the cluster ID, look up the gene symbol
     * of the canonical transcript.  If the symbol is in HUGO, then look up the 
     * Entrez Gene ID */
    geneSymbol = clusterIdToGeneSymbol(conn, clusterId);
    if (isSymbolInHugo(conn, geneSymbol)) {
	hasHugoSymbol = TRUE;
    } else {
	/* If the symbol is not in HUGO, then see if it is the synonym or 
	 * previous symbol of some symbol in HUGO */
	newGeneSymbol = geneSymbolFromSynonymOrPrvSymbol(conn, geneSymbol);
	if (newGeneSymbol != NULL) {
	    freeMem(geneSymbol);
	    geneSymbol = newGeneSymbol;
	    if (isSymbolInHugo(conn, geneSymbol)) {
		hasHugoSymbol = TRUE;
	    }
	}
    } 

    /* If we have a HUGO symbol, then get the Entrez ID.  Otherwise, try
     * Plan C, looking up an Entrez ID for this cluster, and getting
     * the symbol from the HUGO table using the Entrez ID */
    if (!hasHugoSymbol) {
	/* Failing that, try looking up the most frequent Entrez ID for
	 * this cluster */
	entrezGeneId = clusterIdToEntrezId(conn, clusterId);
	newGeneSymbol = entrezIdToGeneSymbol(conn, entrezGeneId);
	if (newGeneSymbol != NULL) {
	    freeMem(geneSymbol);		
	    geneSymbol = newGeneSymbol;
	    hasHugoSymbol = TRUE;
	}
    }
    /* Plan D: try to find the symbol using the knownToEnsembl table 
     * and the ensId field in hgnc */
    if (!hasHugoSymbol) {
	newGeneSymbol = clusterIdToGeneSymbolViaEnsembl(conn, knownIsoformTable,
							clusterId);
	if (newGeneSymbol != NULL) {
	    freeMem(geneSymbol);
	    geneSymbol = newGeneSymbol;
	    hasHugoSymbol = TRUE;
	}
    }

    /* Now that we've tried the many ways to possibly get a gene symbol,
     * work on the Entrez ID. Plan A: Via the HUGO table.  
     * Plan B: via the knownToLocusLink table */
    if (entrezGeneId == NULL) {
	entrezGeneId = geneSymbolToEntrezId(conn, geneSymbol);
    }
    if (entrezGeneId == NULL) {
	entrezGeneId = clusterIdToEntrezId(conn, clusterId);
    }

    /* Create the final gene name, which is a concatenation of 
     * whatever gene symbol we have and the entrez ID if we have one. */
    if (entrezGeneId == NULL) {
	symbolBuf = dyStringCreate("%s|?", geneSymbol);
    } else {
	symbolBuf = dyStringCreate("%s|%s", geneSymbol, entrezGeneId);
	freeMem(entrezGeneId);
    }
    freeMem(geneSymbol);
    geneName = dyStringCannibalize(&symbolBuf);
    return(geneName);
}




int scoreThisTranscript(struct sqlConnection *conn, char *transcriptName, char *geneName)
{
    struct dyString *query = dyStringNew(256);
    struct sqlResult *sr;
    char **row;
    int score;

    dyStringPrintf(query, "select r,g,b from kgColor where kgID = '%s'", transcriptName);
    sr = sqlMustGetResult(conn, query->string);
    row = sqlNextRow(sr);
    if ((atoi(row[0]) == 0 || atoi(row[0]) == 12) && (strcmp(row[1], row[0]) == 0) 
	&& (strcmp(row[2], row[0]) == 0)) {
	/* The transcript is based on a curated transcript */
	score = 1000;
    } else if (atoi(row[0]) == 80 && atoi(row[1]) == 80 && atoi(row[2]) == 160) {
	/* The transcript is based on some other refseq */
	score = 700;
    } else if (stringIn("?", geneName) == NULL) {
	/* The gene maps to a HUGO symbol and ENTREZ ID */
	score = 400;
    } else {
	score = 100;
    }
    sqlFreeResult(&sr);
    dyStringFree(&query);
    return(score);
}

void gpBedOutput(struct genePred *gp, FILE *f)
/* Write out part of gp as bed12. */
{
    /* Figure out # of blocks and min/max of area inside start/end */
    int blockCount = 0;
    int newStart = gp->txEnd, newEnd = gp->txStart;
    int size = 0;
    int i;
    for (i=0; i<gp->exonCount; ++i) {
	int exonStart = gp->exonStarts[i];
	int exonEnd = gp->exonEnds[i];
	if (exonStart < exonEnd) {
	    ++blockCount;
	    newStart = min(exonStart, newStart);
	    newEnd = max(exonEnd, newEnd);
	    size += exonEnd - exonStart;
	}
    }

    /* Output first 10 fields of bed12. */
    fprintf(f, "%s\t%d\t%d\t", gp->chrom, newStart, newEnd);
    fprintf(f, "%s\t", gp->name);
    fprintf(f, "%d\t%s\t%d\t%d\t0\t%d\t", gp->score, gp->strand, 
	    gp->cdsStart, gp->cdsEnd, blockCount);

    /* Output blockSizes field */
    for (i=0; i<gp->exonCount; ++i) {
	int exonStart = gp->exonStarts[i];
	int exonEnd = gp->exonEnds[i];
	if (exonStart < exonEnd)
	    fprintf(f, "%d,", exonEnd - exonStart);
    }
    fprintf(f, "\t");

    /* Output chromStarts field */
    for (i=0; i<gp->exonCount; ++i) {
	int exonStart = gp->exonStarts[i];
	int exonEnd = gp->exonEnds[i];
	if (exonStart < exonEnd)
	    fprintf(f, "%d,", exonStart - newStart);
    }
    fprintf(f, "\n");
}



struct range *storeNewExon(int start, int end)
{
    struct range *exon;

    AllocVar(exon);
    exon->start = start;
    exon->end = end;
    exon->val = NULL;
    return(exon);
}

void addComponentExon(int exonStart, int exonEnd, struct rbTree *componentExons)
{
    struct range *overlappingExons, *newExon, *oldExon;
    struct slRef *newExonList, *nextNewExon;

    int oldExonStart, oldExonEnd, overlapStart, overlapEnd, prefixStart, suffixEnd;

    verbose(2, "Adding new component exon with range %d to %d\n", exonStart, exonEnd);
    overlappingExons = rangeTreeAllOverlapping(componentExons, exonStart, exonEnd);
    if (overlappingExons == NULL) {
	/* No overlapping exons, just need to add the new one */
	newExon = rangeTreeAdd(componentExons, exonStart, exonEnd);
	verbose(2, "New component exon with ranges %d to %d, no overlaps\n", 
		newExon->start, newExon->end);
    } else {
	/* This exon overlaps exons that are already in the tree.  If there is
	 * just one overlapping exon, and it has the same coordinates as this
	 * exon, then do nothing. */
	if (overlappingExons->start == exonStart && overlappingExons->end == exonEnd
	    && overlappingExons->next == NULL) {
	    verbose(2,"Already in the tree\n");
	} else {
	    /* This exon overlaps different exons that are already in the tree.
	     * For each overlapping exon, remove it from the tree, and divide its
	     * range into three sub-ranges: before the new exon, overlapping the
	     * new exon, and after the new exon.  Queue these new exons to be added.
	     * This might seem roundabout, and it might seem more straightforward to
	     * add the new exons directly, but it turns out that it would mess up
	     * the internal data state of a non-reentrant function if we were to do so */
	    newExonList = NULL;
	    for (oldExon = overlappingExons; oldExon != NULL; oldExon = oldExon->next) {
		oldExonStart = oldExon->start;
		oldExonEnd = oldExon->end;
		rbTreeRemove(componentExons, oldExon);
		verbose(2, "Removing an exon from %d to %d\n", oldExonStart, oldExonEnd);
		prefixStart = min(exonStart, oldExonStart);
		overlapStart = max(exonStart, oldExonStart);
		overlapEnd = min(exonEnd, oldExonEnd);
		suffixEnd = max(exonEnd, oldExonEnd);
		verbose(2, "Dividing range into: %d to %d, %d to %d, and %d to %d\n",
			prefixStart, overlapStart, overlapStart, overlapEnd, 
			overlapEnd, suffixEnd);
		if (prefixStart < overlapStart)
		    slAddHead(&newExonList, storeNewExon(prefixStart, overlapStart));
		if (overlapStart < overlapEnd)
		    slAddHead(&newExonList, storeNewExon(overlapStart, overlapEnd));
		if (overlapEnd < suffixEnd)
		    slAddHead(&newExonList, storeNewExon(overlapEnd, suffixEnd));
	    }
	    for (nextNewExon = newExonList; nextNewExon != NULL; 
		 nextNewExon = nextNewExon->next) {
		newExon = (struct range *) nextNewExon;
		addComponentExon(newExon->start, newExon->end, componentExons);
	    }
	    slFreeList(&newExonList);
	}
    }
}

void addExons(struct genePred *gp, struct rbTree *componentExons, struct rbTree*compositeExons)
/* Add the exons from the transcripts to the range trees */
{
    int ii;
    struct range *newExon;

    for (ii = 0; ii < gp->exonCount; ii++) {
	/* Start with a new composite exon with the bounds of exon[ii].  If there are any
	 * overlapping composite exons, update the bounds of this exon to reflect the maximal
	 * span of the overlap, and remove the overlapping exon from the tree.  When done, add
	 * this new exon to the range tree.  Most of this is done implicitly by rangeTreeAdd. */
	newExon = rangeTreeAdd(compositeExons, gp->exonStarts[ii], gp->exonEnds[ii]);
	verbose(2, "New composite exon on %s with ranges %d to %d, final ranges %d to %d\n", gp->chrom, 
	  gp->exonStarts[ii], gp->exonEnds[ii], newExon->start, newExon->end); 
	/* Start with a new component exon with the bounds of exon[ii].  If there are any
	 * intersecting exons, remove them and replace them with three new exons: prefix, overlap,
	 * and suffix.  Do this recursively, reflecting that each new exon added (prefix, overlap,
	 * suffix) might itself intersect some other exon */
	addComponentExon(gp->exonStarts[ii], gp->exonEnds[ii], componentExons);
    }
}

void outputGeneModel(int clusterId, struct rbTree *compositeExons, char *chrom,
		     char strand, char *geneName, int maxGeneScore,
		     FILE *geneModelFile)
{
    struct slRef *exonList, *thisExon;
    struct range *exonBounds;
    struct bed *bed;
    int ii, exonCount;
    char buffer[255];

    AllocVar(bed);
    exonList = rbTreeItems(compositeExons);
    exonCount = 0;
    for (thisExon = exonList; thisExon != NULL; thisExon = thisExon->next) {
	exonCount++;
    }
    bed->next = NULL;
    bed->chrom = cloneString(chrom);
    safef(buffer, sizeof(buffer), "%s;%d", geneName, clusterId);
    bed->name = cloneString(buffer);
    bed->score = maxGeneScore;
    bed->strand[0] = strand;
    bed->strand[1] = '\0';
    bed->itemRgb = 0;
    bed->blockCount = exonCount;
    bed->blockSizes = needMem(exonCount * sizeof(int));
    bed->chromStarts = needMem(exonCount * sizeof(int));
    for (thisExon = exonList, ii = 0; thisExon != NULL; thisExon = thisExon->next, ii++) {
	exonBounds = thisExon->val;
	if (ii == 0) {
	    bed->chromStart = bed->thickStart = exonBounds->start;
	}
	bed->blockSizes[ii] = exonBounds->end - exonBounds->start;
	bed->chromStarts[ii] = exonBounds->start - bed->chromStart;
    }
    bed->chromEnd = bed->thickEnd = exonBounds->end;
    bedTabOutN(bed, 12, geneModelFile);
    bedFree(&bed);
    slFreeList(&exonList);
}



struct bed *newExonBed(struct range *exonBounds, char *chrom, char strand, 
		       int clusterId)
{
    struct bed *bed;
    struct dyString *dyTmp;

    AllocVar(bed);
    bed->next = NULL;
    bed->chrom = cloneString(chrom);
    bed->chromStart = exonBounds->start;
    bed->chromEnd = exonBounds->end;
    dyTmp = dyStringCreate("%s:%d-%d:%c;%d", chrom, exonBounds->start + 1,
			   exonBounds->end, strand, clusterId);
    bed->name = dyStringCannibalize(&dyTmp);
    bed->score = 1000;
    bed->strand[0] = strand;
    bed->strand[1] = '\0';
    makeItBed12(bed, 6);
    return(bed);
}


void outputExonSet(struct rbTree *compositeExons, char *chrom,
		   char strand, int clusterId, FILE *exonFile)
{
    struct slRef *exonList, *thisExon;
    struct range *exonBounds;
    struct bed *thisExonBed;
    int exonId = 1;

    exonList = rbTreeItems(compositeExons);
    for (thisExon = exonList; thisExon != NULL; thisExon = thisExon->next) {
	exonBounds = thisExon->val;
	thisExonBed = newExonBed(exonBounds, chrom, strand, clusterId);
	bedTabOutN(thisExonBed, 12, exonFile);
	bedFree(&thisExonBed);
	exonId++;
    }
    slFreeList(&exonList);
}


void oneChromStrand(struct sqlConnection *conn1, struct sqlConnection *conn2, 
		    char *transcripts, char *transcriptClusters, 
		    FILE *componentExonFile, FILE *compositeExonFile, 
		    FILE *transcriptFile, FILE *geneModelFile,
		    char *chromName, char strand)
{
    struct sqlResult *sr;
    char **row;
    struct hash *geneStructHash;
    struct genePred *gp;
    struct dyString *query = dyStringNew(256);
    int geneClusterId, prevGeneClusterId;
    char *geneName, *transcriptName;
    struct rbTree *componentExons = NULL;
    struct rbTree *compositeExons = NULL;
    int maxGeneScore;
    char buffer[255];

    /* Get the list of transcript clusters for this chrom and strand.  Each transcript
     * cluster entry contains the cluster ID and transcript name.  Transcripts with the 
     * same cluster ID are part of the same gene. */
    geneStructHash = transcriptsThisChromAndStrand(conn1, transcripts, chromName, strand);
    prevGeneClusterId = -1;
    dyStringPrintf(query, "select ki.* from %s ki, %s kg where ki.transcript = kg.name and kg.chrom = '%s' and kg.strand = '%c' order by ki.clusterId", 
		   transcriptClusters, transcripts, chromName, strand);
    sr = sqlGetResult(conn1, query->string);
    while ((row = sqlNextRow(sr)) != NULL) {
	geneClusterId = atoi(row[0]);
	if (geneClusterId != prevGeneClusterId) {
	    /* Starting on a new gene, output data from the old one */
	    if (prevGeneClusterId > 0) {
		outputGeneModel(prevGeneClusterId, compositeExons, chromName, strand, 
				geneName, maxGeneScore, geneModelFile);
		outputExonSet(compositeExons, chromName, strand, prevGeneClusterId,
			      compositeExonFile);
		rangeTreeFree(&compositeExons);
		outputExonSet(componentExons, chromName, strand, prevGeneClusterId,
			      componentExonFile);
		rangeTreeFree(&componentExons);
	    }
	    prevGeneClusterId = geneClusterId;
	    geneName = geneNameThisCluster(conn2, geneClusterId, transcriptClusters);
	    maxGeneScore = 0;
	    componentExons = rangeTreeNew();
	    compositeExons = rangeTreeNew();
	} 
	transcriptName = row[1];
	gp = hashMustFindVal(geneStructHash, transcriptName);
	gp->score = scoreThisTranscript(conn2, transcriptName, geneName);
	safef(buffer, sizeof(buffer), "%s;%d", transcriptName, geneClusterId);
	gp->name = cloneString(buffer);
	maxGeneScore = max(gp->score, maxGeneScore);
	gpBedOutput(gp, transcriptFile);
	addExons(gp, componentExons, compositeExons);
    }
    if (compositeExons != NULL) {
	outputGeneModel(geneClusterId, compositeExons, chromName, strand, 
			geneName, maxGeneScore, geneModelFile);
	outputExonSet(compositeExons, chromName, strand, geneClusterId, 
		      compositeExonFile);
	rangeTreeFree(&compositeExons);
    }
    if (componentExons != NULL) {
	outputExonSet(componentExons, chromName, strand, geneClusterId, 
		      componentExonFile);
	rangeTreeFree(&componentExons);
    }
    dyStringFree(&query);
    sqlFreeResult(&sr);
    hashFree(&geneStructHash);
}

void transcriptsToGenesAndExons(char *database, char *transcripts, char *transcriptClusters, 
				FILE *componentExonFile, FILE *compositeExonFile, 
				FILE *transcriptFile, FILE *geneModelFile)
{
    struct slName *chromList, *chrom;
    struct sqlConnection *conn1 = sqlConnect(database); 
    struct sqlConnection *conn2 = sqlConnect(database); 
    chromList = hAllChromNames(database);
    for (chrom = chromList; chrom != NULL; chrom = chrom->next) {
	oneChromStrand(conn1, conn2, transcripts, transcriptClusters, 
		       componentExonFile, compositeExonFile, transcriptFile, 
		       geneModelFile, chrom->name, '+');
	oneChromStrand(conn1, conn2, transcripts, transcriptClusters, 
		       componentExonFile, compositeExonFile, transcriptFile, 
		       geneModelFile, chrom->name, '-');  
    }
}


int main(int argc, char *argv[])
{
    FILE *componentExonFile;
    FILE *compositeExonFile;
    FILE *transcriptFile;
    FILE *geneModelFile;

    optionInit(&argc, argv, options);
    if (argc != 8)
	usage();
    componentExonFile = mustOpen(argv[4], "w");
    compositeExonFile = mustOpen(argv[5], "w");
    transcriptFile = mustOpen(argv[6], "w");
    geneModelFile = mustOpen(argv[7], "w");
    transcriptsToGenesAndExons(argv[1], argv[2], argv[3], componentExonFile,
			       compositeExonFile, transcriptFile, geneModelFile);
    fclose(componentExonFile);
    fclose(compositeExonFile);
    fclose(transcriptFile);
    fclose(geneModelFile);
    exit(0);
}
