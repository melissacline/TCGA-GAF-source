
gaf21File = /hive/users/cline/TCGA/GAF2.1/TCGA.hg19.June2011.with_dbSNP.gaf.gz
inputDir = data/hg19-init
scratchDir = data/scratch
gafDir = data/gaf
outputBedDir = data/hg19-final
testDir = data/test
testInput = ${testDir}/input
testOutput = ${testDir}/output
transcriptTestFields = 2-5,8-9,11-13
geneTestFields = 2-5,7,9-10,12-17


geneSetContents = ${gafDir}/gene.genome.gaf ${gafDir}/transcript.genome.gaf \
	${gafDir}/compositeExon.genome.gaf ${gafDir}/componentExon.genome.gaf \
	${gafDir}/junction.genome.gaf \
	${gafDir}/compositeExon.gene.gaf ${gafDir}/compositeExon.transcript.gaf \
	${gafDir}/componentExon.gene.gaf ${gafDir}/componentExon.transcript.gaf \
	${gafDir}/junction.transcript.gaf

miRnaSetContents = ${gafDir}/pre-miRNA.genome.gaf ${gafDir}/miRNA.genome.gaf \
	${gafDir}/miRNA.pre-miRNA.gaf

all:	${gafDir}/geneSet.gaf \
	${outputBedDir}/gene.genome.bb ${outputBedDir}/compositeExon.genome.bb \
	${outputBedDir}/componentExon.genome.bb ${outputBedDir}/transcript.genome.bb \
	${outputBedDir}/junction.genome.bb \
	${testOutput}/transcript.genome.diff ${testOutput}/gene.genome.diff \
	${testOutput}/transcript.gene.diff ${gafDir}/transcript.gene.gaf \
	${testOutput}/compositeExon.transcript.diff \
	${testOutput}/componentExon.transcript.diff \
	${testOutput}/junction.transcript.diff \
	${gafDir}/miRnaSet.gaf \
	${outputBedDir}/pre-miRNA.genome.bb ${outputBedDir}/miRNA.genome.bb \
        ${testOutput}/pre-miRNA.genome.diff ${testOutput}/miRNA.genome.diff \
        ${testOutput}/miRNA.pre-miRNA.diff \
	${outputBedDir}/affySnp.genome.bb ${outputBedDir}/dbSNP.genome.bb \
	${outputBedDir}/dbSNP.genome.bb ${testOutput}/dbSNP.genome.diff \
        ${outputBedDir}/MAprobe.genome.bb ${testOutput}/MAprobe.genome.diff \
        ${outputBedDir}/AffySNP.genome.bb

${gafDir}/geneSet.gaf:	${geneSetContents}
	cat ${geneSetContents} | awk '{$$1 = NR; print}' > $@

${testOutput}/gene.genome.diff:	${testInput}/gene.genome.2.1.gaf ${testInput}/gene.genome.3.0.gaf
	cat ${testInput}/gene.genome.2.1.gaf \
	| awk -F'\t' '{ print $$2, $$3, $$4, $$5, $$6, $$7, $$9, $$10, $$12, $$13, $$14, $$15, $$16, $$17 }' \
        > ${scratchDir}/gene.genome.2.1.subset
	cat ${testInput}/gene.genome.3.0.gaf \
	| awk -F'\t' '{ print $$2, $$3, $$4, $$5, $$6, $$7, $$9, $$10, $$12, $$13, $$14, $$15, $$16, $$17 }' \
        > ${scratchDir}/gene.genome.3.0.subset
	diff ${scratchDir}/gene.genome.2.1.subset ${scratchDir}/gene.genome.3.0.subset > $@

${testInput}/gene.genome.3.0.gaf:	${gafDir}/gene.genome.gaf ${testDir}/testGenes.txt
	cat ${testDir}/testGenes.txt \
	| awk '{ print "grep \"" $$1 "|" $$2 "\" ${gafDir}/gene.genome.gaf"}' |bash \
	> $@

${testInput}/gene.genome.2.1.gaf: ${scratchDir}/gene.genome.gaf21.gaf
	cat ${testDir}/testGenes.txt \
	| awk '{ print "grep \"" $$1 "|" $$2  "\" ${scratchDir}/gene.genome.gaf21.gaf"}' |bash > $@

${scratchDir}/gene.genome.gaf21.gaf: ${testDir}/testGenes.txt
	zcat ${gaf21File} \
	| awk -F'\t' '$$3 == "gene" && $$9 == "genome" { print }' \
	> ${scratchDir}/gene.genome.gaf21.gaf 

${gafDir}/gene.genome.gaf:	${inputDir}/gene.genome.bed
	liftOver $< data/GRCh37-lite/hg19.GRCh37-lite.over.chain ${scratchDir}/gene.genome.preGaf.GRCh37-lite.bed /dev/null
	hgsql hg19 < sql/gafGeneXref.sql
	scripts/geneToGeneXref.py $< ${scratchDir}/gene.genome.preGaf.GRCh37-lite.bed
	scripts/makeGeneSet.py ${scratchDir}/gene.genome.preGaf.GRCh37-lite.bed |sort -k2,2 > $@ 


${testOutput}/transcript.genome.diff:	${testInput}/transcript.genome.2.1.gaf ${testInput}/transcript.genome.3.0.gaf
	cat ${testInput}/transcript.genome.2.1.gaf \
	| awk -F'\t' '{ print $$2, $$3, $$4, $$5, $$9, $$10, $$12, $$13, $$14, $$15 }' \
        > ${scratchDir}/transcript.genome.2.1.subset
	cat ${testInput}/transcript.genome.3.0.gaf \
	| awk -F'\t' '{ print $$2, $$3, $$4, $$5, $$9, $$10, $$12, $$13, $$14, $$15 }' \
        > ${scratchDir}/transcript.genome.3.0.subset
	diff ${scratchDir}/transcript.genome.2.1.subset ${scratchDir}/transcript.genome.3.0.subset > $@


${testInput}/transcript.genome.3.0.gaf:	${gafDir}/transcript.genome.gaf ${testDir}/testTranscripts.txt
	cat ${testDir}/testTranscripts.txt \
	| awk '{ print "grep", $$1, "${gafDir}/transcript.genome.gaf"}' |bash > $@


${testInput}/transcript.genome.2.1.gaf:	${testDir}/testTranscripts.txt ${scratchDir}/transcript.genome.gaf21.gaf
	cat ${testDir}/testTranscripts.txt \
	| awk '{ print "grep", $$1, "${scratchDir}/transcript.genome.gaf21.gaf"}'  | bash > $@

${scratchDir}/transcript.genome.gaf21.gaf:
	zcat ${gaf21File} \
	| awk -F'\t' '$$3 == "transcript" && $$9 == "genome" { print }' > $@


${gafDir}/transcript.genome.gaf:	${inputDir}/transcript.genome.bed ${gafDir}/gene.genome.gaf
	liftOver $< data/GRCh37-lite/hg19.GRCh37-lite.over.chain ${scratchDir}/transcript.genome.preGaf.GRCh37-lite.bed /dev/null
	scripts/makeTranscriptSet.py ${scratchDir}/transcript.genome.preGaf.GRCh37-lite.bed |sort > $@ 


${testOutput}/transcript.gene.diff:	${testInput}/transcript.gene.2.1.gaf ${testInput}/transcript.gene.3.0.gaf
	cat ${testInput}/transcript.gene.2.1.gaf \
	| awk -F'\t' '{ print $$2, $$3, $$4, $$5, $$8, $$9, $$10, $$11, $$12, $$13, $$14, $$15, $$16, $$17, $$18 }' \
        > ${testInput}/transcript.gene.2.1.subset
	cat ${testInput}/transcript.gene.3.0.gaf \
	| awk -F'\t' '{ print $$2, $$3, $$4, $$5, $$8, $$9, $$10, $$11, $$12, $$13, $$14, $$15, $$16, $$17, $$18 }' \
        > ${testInput}/transcript.gene.3.0.subset
	diff ${testInput}/transcript.gene.2.1.subset ${testInput}/transcript.gene.3.0.subset > $@

${testInput}/transcript.gene.3.0.gaf:	${testInput}/transcript.genome.3.0.gaf ${testInput}/gene.genome.3.0.gaf
	scripts/featureToComposite.py ${testInput}/transcript.genome.3.0.gaf ${testInput}/gene.genome.3.0.gaf |sort -k2,2 > $@

${testInput}/transcript.gene.2.1.gaf: ${scratchDir}/transcript.gene.2.1.gaf ${testDir}/testTranscripts.txt
	cat ${testDir}/testTranscripts.txt |awk '{print "grep", $$1, "${scratchDir}/transcript.gene.2.1.gaf"}' |bash > ${scratchDir}/transcript.gene.2.1.superset.gaf
	cat ${testDir}/testGenes.txt |awk '{ print "grep -w \"" $$1 "|" $$2  "\" ${scratchDir}/transcript.gene.2.1.superset.gaf"}' |bash |sort -k2,2 > $@
	rm ${scratchDir}/transcript.gene.2.1.superset.gaf

${scratchDir}/transcript.gene.2.1.gaf:
	zcat ${gaf21File} \
	| awk -F'\t' '$$3 == "transcript" && $$9 == "gene" { print }' > $@

${gafDir}/transcript.gene.gaf: ${gafDir}/transcript.genome.gaf ${gafDir}/gene.genome.gaf
	scripts/featureToComposite.py ${gafDir}/transcript.genome.gaf ${gafDir}/gene.genome.gaf > $@

${testOutput}/componentExon.transcript.diff:	${testInput}/componentExon.transcript.2.1.gaf ${testInput}/componentExon.transcript.3.0.gaf
	cat ${testInput}/componentExon.transcript.2.1.gaf \
	| awk -F'\t' '{ print $$2, $$3, $$4, $$5, $$8, $$9, $$10, $$11, $$13, $$14, $$15, $$16, $$17, $$18 }' \
        > ${testInput}/componentExon.transcript.2.1.subset
	cat ${testInput}/componentExon.transcript.3.0.gaf \
	| awk -F'\t' '{ print $$2, $$3, $$4, $$5, $$8, $$9, $$10, $$11, $$13, $$14, $$15, $$16, $$17, $$18 }' \
        > ${testInput}/componentExon.transcript.3.0.subset
	diff ${testInput}/componentExon.transcript.2.1.subset ${testInput}/componentExon.transcript.3.0.subset > $@

${testInput}/componentExon.transcript.3.0.gaf:	${testInput}/componentExon.genome.3.0.gaf ${testInput}/transcript.genome.3.0.gaf
	scripts/featureToComposite.py ${testInput}/componentExon.genome.3.0.gaf ${testInput}/transcript.genome.3.0.gaf |sort -k2,2 |uniq > $@

${testInput}/componentExon.genome.3.0.gaf:	${gafDir}/componentExon.genome.gaf ${testDir}/testGenes.txt
	cat ${testDir}/testGenes.txt \
	| awk '{ print "grep", $$1, "${gafDir}/componentExon.genome.gaf"}' |bash |sort -k2,2 |uniq > $@

${testInput}/componentExon.transcript.2.1.gaf: ${scratchDir}/componentExon.transcript.2.1.gaf ${testDir}/testTranscripts.txt ${testDir}/testGenes.txt
	cat ${testDir}/testTranscripts.txt \
	|awk '{print "grep", $$1, "${scratchDir}/componentExon.transcript.2.1.gaf"}'  \
	|bash > ${scratchDir}/componentExon.transcript.2.1.superset.gaf
	cat ${testDir}/testGenes.txt |awk '{ print "grep -w \"" $$1 "|" $$2  "\" ${scratchDir}/componentExon.transcript.2.1.superset.gaf"}' \
	|bash |sort -k2,2 |uniq > $@
	rm ${scratchDir}/componentExon.transcript.2.1.superset.gaf



${scratchDir}/componentExon.transcript.2.1.gaf:
	zcat ${gaf21File} \
	| awk -F'\t' '$$3 == "componentExon" && $$9 == "transcript" { print }' > $@

${gafDir}/componentExon.transcript.gaf: ${gafDir}/componentExon.genome.gaf ${gafDir}/transcript.genome.gaf
	scripts/featureToComposite.py ${gafDir}/componentExon.genome.gaf ${gafDir}/transcript.genome.gaf > $@

${gafDir}/componentExon.gene.gaf: ${gafDir}/componentExon.genome.gaf ${gafDir}/gene.genome.gaf
	scripts/featureToComposite.py ${gafDir}/componentExon.genome.gaf ${gafDir}/gene.genome.gaf > $@

${testOutput}/compositeExon.transcript.diff:	${testInput}/compositeExon.transcript.2.1.gaf ${testInput}/compositeExon.transcript.3.0.gaf
	cat ${testInput}/compositeExon.transcript.2.1.gaf \
	| awk -F'\t' '{ print $$2, $$3, $$4, $$5, $$8, $$9, $$10, $$11, $$13, $$14, $$15, $$16, $$17, $$18, $$19 }' \
        > ${testInput}/compositeExon.transcript.2.1.subset
	cat ${testInput}/compositeExon.transcript.3.0.gaf \
	| awk -F'\t' '{ print $$2, $$3, $$4, $$5, $$8, $$9, $$10, $$11, $$13, $$14, $$15, $$16, $$17, $$18, $$19 }' \
        > ${testInput}/compositeExon.transcript.3.0.subset
	diff ${testInput}/compositeExon.transcript.2.1.subset ${testInput}/compositeExon.transcript.3.0.subset > $@

${testInput}/compositeExon.transcript.3.0.gaf:	${testInput}/compositeExon.genome.3.0.gaf ${testInput}/transcript.genome.3.0.gaf
	scripts/featureToComposite.py ${testInput}/compositeExon.genome.3.0.gaf ${testInput}/transcript.genome.3.0.gaf |sort -k2,2 |uniq > $@

${gafDir}/compositeExon.gene.gaf: ${gafDir}/compositeExon.genome.gaf ${gafDir}/gene.genome.gaf
	scripts/featureToComposite.py ${gafDir}/compositeExon.genome.gaf ${gafDir}/gene.genome.gaf > $@

${testInput}/compositeExon.genome.3.0.gaf:	${gafDir}/compositeExon.genome.gaf ${testDir}/testGenes.txt
	cat ${testDir}/testGenes.txt \
	| awk '{ print "grep", $$1, "${gafDir}/compositeExon.genome.gaf"}' \
        | bash |sort -k2,2 |uniq > $@


${testInput}/compositeExon.transcript.2.1.gaf: ${scratchDir}/compositeExon.transcript.2.1.gaf ${testDir}/testTranscripts.txt ${testDir}/testGenes.txt
	cat ${testDir}/testTranscripts.txt \
	|awk '{print "grep", $$1, "${scratchDir}/compositeExon.transcript.2.1.gaf"}'  \
	|bash > ${scratchDir}/compositeExon.transcript.2.1.superset.gaf
	cat ${testDir}/testGenes.txt |awk '{ print "grep -w \"" $$1 "|" $$2  "\" ${scratchDir}/compositeExon.transcript.2.1.superset.gaf"}' \
	|bash |sort -k2,2 > $@
	rm ${scratchDir}/compositeExon.transcript.2.1.superset.gaf

${scratchDir}/compositeExon.transcript.2.1.gaf:
	zcat ${gaf21File} \
	| awk -F'\t' '$$3 == "compositeExon" && $$9 == "transcript" { print }' > $@

${gafDir}/compositeExon.transcript.gaf: ${gafDir}/compositeExon.genome.gaf ${gafDir}/transcript.genome.gaf
	scripts/featureToComposite.py ${gafDir}/compositeExon.genome.gaf ${gafDir}/transcript.genome.gaf > $@

${testOutput}/junction.transcript.diff:	${testInput}/junction.transcript.2.1.gaf ${testInput}/junction.transcript.3.0.gaf
	cat ${testInput}/junction.transcript.2.1.gaf \
	| awk -F'\t' '{ print $$2, $$3, $$4, $$5, $$8, $$9, $$10, $$11, $$13, $$14, $$15, $$16, $$17, $$18, $$19 }' \
        |sort |uniq > ${testInput}/junction.transcript.2.1.subset
	cat ${testInput}/junction.transcript.3.0.gaf \
	| awk -F'\t' '{ print $$2, $$3, $$4, $$5, $$8, $$9, $$10, $$11, $$13, $$14, $$15, $$16, $$17, $$18, $$19 }' \
        |sort |uniq > ${testInput}/junction.transcript.3.0.subset
	diff ${testInput}/junction.transcript.2.1.subset ${testInput}/junction.transcript.3.0.subset > $@

${testInput}/junction.transcript.3.0.gaf:	${testInput}/junction.genome.3.0.gaf ${testInput}/transcript.genome.3.0.gaf 
	cat ${testDir}/testTranscripts.txt \
	|awk '{print "grep", $$1, "${gafDir}/junction.transcript.gaf"}'  \
	|bash > ${scratchDir}/junction.transcript.3.0.superset.gaf
	cat ${testDir}/testGenes.txt |awk '{ print "grep -w \"" $$1 "|" $$2  "\" ${scratchDir}/junction.transcript.3.0.superset.gaf"}' \
	|bash |sort -k2,2 |uniq> $@
	rm ${scratchDir}/junction.transcript.3.0.superset.gaf

${testInput}/junction.transcript.2.1.gaf: ${scratchDir}/junction.transcript.2.1.gaf ${testDir}/testTranscripts.txt
	cat ${testDir}/testTranscripts.txt \
	|awk '{print "grep", $$1, "${scratchDir}/junction.transcript.2.1.gaf"}'  \
	|bash > ${scratchDir}/junction.transcript.2.1.superset.gaf
	cat ${testDir}/testGenes.txt |awk '{ print "grep -w \"" $$1 "|" $$2  "\" ${scratchDir}/junction.transcript.2.1.superset.gaf"}' \
	|bash |sort -k2,2 |uniq > $@
	rm ${scratchDir}/junction.transcript.2.1.superset.gaf

${scratchDir}/junction.transcript.2.1.gaf:
	zcat ${gaf21File} \
	| awk -F'\t' '$$3 == "junction" && $$9 == "transcript" { print }' > $@

${gafDir}/junction.transcript.gaf: ${gafDir}/junction.genome.gaf ${gafDir}/transcript.genome.gaf
	scripts/junctionToTranscript.py ${gafDir}/junction.genome.gaf ${gafDir}/transcript.genome.gaf > $@


${testDir}/testGenes.txt ${testDir}/testTranscripts.txt: sql/geneTestSet.sql
	rm ${testDir}/testGenes.txt ${testDir}/testTranscripts.txt
	hgsql hg19 < sql/geneTestSet.sql

${outputBedDir}/%.bb:	${gafDir}/%.gaf
	scripts/gafToBed.py $< |sort -k1,1 -k2,2n > ${scratchDir}/$*.postGaf.GRCh37-lite.bed
	liftOver ${scratchDir}/$*.postGaf.GRCh37-lite.bed data/GRCh37-lite/GRCh37-lite.hg19.over.chain ${scratchDir}/$*.postGaf.hg19.bed /dev/null
	bedToBigBed ${scratchDir}/$*.postGaf.hg19.bed /hive/data/genomes/hg19/chrom.sizes $@

${gafDir}/componentExon.genome.gaf:	${inputDir}/componentExon.genome.bed ${gafDir}/gene.genome.gaf
	liftOver $< data/GRCh37-lite/hg19.GRCh37-lite.over.chain ${scratchDir}/componentExon.genome.preGaf.GRCh37-lite.bed /dev/null
	scripts/makeExonSet.py -t componentExon ${scratchDir}/componentExon.genome.preGaf.GRCh37-lite.bed |sort > $@ 

${gafDir}/compositeExon.genome.gaf:	${inputDir}/compositeExon.genome.bed ${gafDir}/gene.genome.gaf
	liftOver $< data/GRCh37-lite/hg19.GRCh37-lite.over.chain ${scratchDir}/compositeExons.preGaf.GRCh37-lite.bed /dev/null
	scripts/makeExonSet.py -t compositeExon ${scratchDir}/compositeExons.preGaf.GRCh37-lite.bed |sort >  $@ 


${gafDir}/junction.genome.gaf:	${inputDir}/transcript.genome.bed ${gafDir}/gene.genome.gaf
	liftOver ${inputDir}/transcript.genome.bed data/GRCh37-lite/hg19.GRCh37-lite.over.chain ${scratchDir}/transcripts.preGaf.GRCh37-lite.bed /dev/null
	scripts/makeJunctionGenome.py  ${scratchDir}/transcripts.preGaf.GRCh37-lite.bed |sort |uniq >  $@ 
	rm ${scratchDir}/transcripts.preGaf.GRCh37-lite.bed

${inputDir}/componentExon.genome.bed ${inputDir}/compositeExon.genome.bed ${inputDir}/transcript.genome.bed ${inputDir}/gene.genome.bed:	src/transcriptsToGenesAndExons
	src/transcriptsToGenesAndExons \
	hg19 knownGene knownIsoformsClustersMerged \
	${inputDir}/componentExon.genome.bed ${inputDir}/compositeExon.genome.bed \
	${inputDir}/transcript.genome.bed ${inputDir}/gene.genome.bed 


#
# miRNA data
#

${gafDir}/miRnaSet.gaf:	${miRnaSetContents}
	cat ${miRnaSetContents} | awk '{$$1 = NR; print}' > $@

#
# When comparing the new pre-miRNA data to the old, don't look at the strand because
# miRBase has changed the strands for many miRNAs.  Also don't look at the gene name
# or gene locus field, because no gene assignment was made in a lot of the old
# miRBase records.
${testOutput}/pre-miRNA.genome.diff:	${testInput}/pre-miRNA.genome.2.1.gaf ${testInput}/pre-miRNA.genome.3.0.gaf 
	cat ${testInput}/pre-miRNA.genome.2.1.gaf \
	| awk -F'\t' '{ print $$2, $$3, $$4, $$9, $$10, $$12, $$13, $$14, $$15}' \
	| sed 's/:+//' |sed 's/:-//' > ${scratchDir}/pre-miRNA.genome.2.1.subset
	cat ${testInput}/pre-miRNA.genome.3.0.gaf \
	| awk -F'\t' '{ print $$2, $$3, $$4, $$9, $$10, $$12, $$13, $$14, $$15}' \
	| sed 's/:+//' |sed 's/:-//' > ${scratchDir}/pre-miRNA.genome.3.0.subset
	diff ${scratchDir}/pre-miRNA.genome.2.1.subset ${scratchDir}/pre-miRNA.genome.3.0.subset > $@

${testInput}/pre-miRNA.genome.2.1.gaf:	${testDir}/testPreMiRna.txt ${scratchDir}/pre-miRNA.genome.gaf21.gaf
	cat ${testDir}/testPreMiRna.txt \
	| awk '{ print "grep \"" $$1 "|\" ${scratchDir}/pre-miRNA.genome.gaf21.gaf"}'  | bash > $@

${testInput}/pre-miRNA.genome.3.0.gaf: ${testDir}/testPreMiRna.txt ${gafDir}/pre-miRNA.genome.gaf
	cat ${testDir}/testPreMiRna.txt \
	| awk '{ print "grep \"" $$1 "|\" ${gafDir}/pre-miRNA.genome.gaf"}'  | bash > $@

${scratchDir}/pre-miRNA.genome.gaf21.gaf:
	zcat ${gaf21File} \
	| awk -F'\t' '$$3 == "pre-miRNA" && $$9 == "genome" { print }' > $@

#
# After making the pre-miRNA GAF data, add the genes to the 
${gafDir}/pre-miRNA.genome.gaf:	${inputDir}/pre-miRNA.genome.bed
	liftOver ${inputDir}/pre-miRNA.genome.bed data/GRCh37-lite/hg19.GRCh37-lite.over.chain ${scratchDir}/pre-miRNA.genome.preGaf.GRCh37-lite.bed /dev/null
	scripts/makePreMiRna.py ${scratchDir}/pre-miRNA.genome.preGaf.GRCh37-lite.bed data/miRNA.dat > $@
	scripts/preMiRnaToGeneXref.py $< $@ 

${inputDir}/pre-miRNA.genome.bed:	data/hsa.gff3
	scripts/gffToBed.py -t miRNA_primary_transcript -n accession_number $< > $@

#
# The list of test pre-miRNAs is created as follows:
#
# awk '{ split($2, tokens, "|"); print tokens[1]}' < ../gaf/pre-miRNA.genome.gaf \
#   | sort > preMiRna.old.txt
# awk '{ split($2, tokens, "|"); print tokens[1]}' \
#   < ../scratch/pre-miRNA.genome.gaf21.gaf |sort > preMiRna.new.txt
# join preMiRna.old.txt preMiRna.new.txt |wc
#    889     889   11492
# join preMiRna.old.txt preMiRna.new.txt |head -250 > testPreMiRna.txt


#
# When looking at the miRNA.genome records, don't look at the miRBase ID 
# (the part that begins 'hsa'), or the strand, because those have been changed
# by miRBase in many cases.  Also don't look at cases where multiple miRNAs
# have been combined into a single entry, because there is nothing to ensure
# that the order of the entries is unchanged since the last GAF build.
${testOutput}/miRNA.genome.diff:	${testInput}/miRNA.genome.2.1.gaf ${testInput}/miRNA.genome.3.0.gaf 
	cat ${testInput}/miRNA.genome.2.1.gaf \
	| awk -F'\t' '{ split($$2, tokens, "|"); print tokens[2], $$3, $$4, $$9, $$10, $$12, $$13, $$14, $$15, $$16, $$17}' \
	|sed 's/:-//g' |sed s'/:+//g' |grep -v ";" |sort \
        > ${scratchDir}/miRNA.genome.2.1.subset
	cat ${testInput}/miRNA.genome.3.0.gaf \
	| awk -F'\t' '{ split($$2, tokens, "|"); print tokens[2], $$3, $$4, $$9, $$10, $$12, $$13, $$14, $$15, $$16, $$17}' \
	|sed 's/:-//g' |sed s'/:+//g' |grep -v ";" |sort \
        > ${scratchDir}/miRNA.genome.3.0.subset
	diff ${scratchDir}/miRNA.genome.2.1.subset ${scratchDir}/miRNA.genome.3.0.subset > $@

${testInput}/miRNA.genome.2.1.gaf:	${testDir}/testPreMiRna.txt ${scratchDir}/miRNA.genome.gaf21.gaf
	cat ${testDir}/testPreMiRna.txt \
	| awk '{ print "grep \"" $$1 "|\" ${scratchDir}/miRNA.genome.gaf21.gaf"}'  | bash > $@

${testInput}/miRNA.genome.3.0.gaf: ${testDir}/testPreMiRna.txt ${gafDir}/miRNA.genome.gaf
	cat ${testDir}/testPreMiRna.txt \
	| awk '{ print "grep \"" $$1 "|\" ${gafDir}/miRNA.genome.gaf"}'  | bash > $@

${scratchDir}/miRNA.genome.gaf21.gaf:
	zcat ${gaf21File} \
	| awk -F'\t' '$$3 == "miRNA" && $$9 == "genome" { print }' > $@


${gafDir}/miRNA.genome.gaf:	${inputDir}/miRNA.genome.bed 
	scripts/makeMiRna.py $< data/hsa.gff3 |sed 's/|/ /' |sort -k3,13\
        |sed 's/ /|/' > ${scratchDir}/miRNA.genome.uncombined.gaf
	scripts/combineMiRnas.py ${scratchDir}/miRNA.genome.uncombined.gaf > $@

${inputDir}/miRNA.genome.bed:	data/hsa.gff3
	scripts/gffToBed.py -t miRNA -n ID $< > $@

#
# When looking at the miRNA.genome records, don't look at the following things,
# that are known to have changed since the last release:
# - the strand.  In this case, it appears in the locus string.
# - As consequence, the composite coordinates are often off.  If a miRNA was
#   at the 5' end of a pre-miRNA and the strand changed but the coordinates
#   stayed the same, this miRNA is now at the 3' end of the pre-miRNA.
# - the miRBase name, the part that begins 'hsa'.  These are recognizable from
#   miRBase 15, but many have a new suffix (*-1, *-2) to distinguish paralogs
${testOutput}/miRNA.pre-miRNA.diff:	${testInput}/miRNA.pre-miRNA.2.1.gaf ${testInput}/miRNA.pre-miRNA.3.0.gaf 
	cat ${testInput}/miRNA.pre-miRNA.2.1.gaf \
	| awk -F'\t' '{ split($$2, tokens, "|"); print tokens[2], $$3, $$4, $$8, $$9, $$10, $$13, $$14, $$16, $$17 }' \
	|sed 's/:-//g' |sed s'/:+//g' |grep -v ";" |sort -k1,1 -k8,8 \
	> ${scratchDir}/miRNA.pre-miRNA.2.1.subset
	cat ${testInput}/miRNA.pre-miRNA.3.0.gaf \
	| awk -F'\t' '{ split($$2, tokens, "|"); print tokens[2], $$3, $$4, $$8, $$9, $$10, $$13, $$14, $$16, $$17}' \
	|sed 's/:-//g' |sed s'/:+//g' |grep -v ";" |sort -k1,1 -k8,8 \
	> ${scratchDir}/miRNA.pre-miRNA.3.0.subset
	diff ${scratchDir}/miRNA.pre-miRNA.2.1.subset ${scratchDir}/miRNA.pre-miRNA.3.0.subset > $@

${testInput}/miRNA.pre-miRNA.2.1.gaf:	${testDir}/testPreMiRna.txt ${scratchDir}/miRNA.pre-miRNA.gaf21.gaf
	cat ${testDir}/testPreMiRna.txt \
	| awk '{ print "grep \"" $$1 "|\" ${scratchDir}/miRNA.pre-miRNA.gaf21.gaf"}'  | bash > $@

${testInput}/miRNA.pre-miRNA.3.0.gaf: ${testDir}/testPreMiRna.txt ${gafDir}/miRNA.pre-miRNA.gaf
	cat ${testDir}/testPreMiRna.txt \
	| awk '{ print "grep \"" $$1 "|\" ${gafDir}/miRNA.pre-miRNA.gaf"}'  | bash > $@

${scratchDir}/miRNA.pre-miRNA.gaf21.gaf:
	zcat ${gaf21File} \
	| awk -F'\t' '$$3 == "miRNA" && $$9 == "pre-miRNA" { print }' > $@


${gafDir}/miRNA.pre-miRNA.gaf:	${gafDir}/miRNA.genome.gaf ${gafDir}/pre-miRNA.genome.gaf
	scripts/miRnasToPreMiRnas.py ${scratchDir}/miRNA.genome.uncombined.gaf ${gafDir}/pre-miRNA.genome.gaf > $@

#
# dbSNP
#

${testOutput}/dbSNP.genome.diff: ${testInput}/dbSNP.genome.2.1.gaf ${testInput}/dbSNP.genome.3.0.gaf
	cat ${testInput}/dbSNP.genome.2.1.gaf \
	| awk -F'\t' '{ print $$2, $$3, $$4, $$6, $$7, $$9, $$10, $$12, $$13, $$14, $$15, $$16, $$18, $$19}' > ${scratchDir}/dbSNP.genome.2.1.subset
	cat ${testInput}/dbSNP.genome.3.0.gaf \
	| awk -F'\t' '{ print $$2, $$3, $$4, $$6, $$7, $$9, $$10, $$12, $$13, $$14, $$15, $$16, $$18, $$19}' > ${scratchDir}/dbSNP.genome.3.0.subset
	diff ${scratchDir}/dbSNP.genome.2.1.subset ${scratchDir}/dbSNP.genome.3.0.subset > $@

${testInput}/dbSNP.genome.3.0.gaf:	${inputDir}/dbSNP-test.bed
	liftOver $< data/GRCh37-lite/hg19.GRCh37-lite.over.chain ${scratchDir}/dbSNP-test.genome.preGaf.GRCh37-lite.bed /dev/null
	scripts/makeDbSnp.py ${scratchDir}/dbSNP-test.genome.preGaf.GRCh37-lite.bed $<\
	|sort -k2,2 |scripts/combineSnps.py > $@ 

#${testInput}/dbSNP.genome.2.1.gaf:	${testDir}/testDbSNP.txt ${scratchDir}/dbSNP.genome.gaf21.gaf
#	cat ${testDir}/testDbSNP.txt \
#	| awk '{ print "grep \"" $$1 "|\" ${scratchDir}/dbSNP.gaf21.gaf"}'  | bash > $@

${scratchDir}/dbSNP.genome.gaf21.gaf:
	zcat ${gaf21File} \
	| awk -F'\t' '$$3 == "dbSNP" && $$9 == "genome" { print }' > $@

${gafDir}/dbSNP.genome.gaf:	${inputDir}/dbSNP.bed ${gafDir}/gene.genome.gaf
	exit  # Clean me up!
	liftOver $< data/GRCh37-lite/hg19.GRCh37-lite.over.chain ${scratchDir}/dbSNP.genome.preGaf.GRCh37-lite.bed /dev/null
	scripts/makeDbSnp.py ${scratchDir}/dbSNP.genome.preGaf.GRCh37-lite.bed $< > ${scratchDir}/dbSNP.raw.gaf
	rm ${scratchDir}/dbSNP.genome.preGaf.GRCh37-lite.bed
	cat ${scratchDir}/dbSNP.raw.gaf |sort -k2,2 | scripts/combineSnps.py > $@ 

#
# MAProbe
data/MAprobe.hg19.bed:	${gafDir}/MAprobe.genome.gaf
	scripts/gafToBed.py $< > ${scratchDir}/MAprobe.postGaf.GRCh37-lite.bed
	liftOver ${scratchDir}/MAprobe.postGaf.GRCh37-lite.bed data/GRCh37-lite/GRCh37-lite.hg19.over.chain $@ /dev/null
	hgLoadBed hg19 gafMaProbe $@

#
# When comparing old and new versions of the MAprobe genomic alignments, don't look at:
# - the gene and gene locus fields, because a large proportion have changed.
# - minus strand entries, because there is a bug in the old data such that some
#   entries are given with backwards minus strand coordinates (i.e. coordinates not
#   in increasing order w.r.t. the positive strand)
# - entries with multiple alignments, because there's no guarantee that we'll combine
#   the composite strings in the same order as the old data.
#
${testOutput}/MAprobe.genome.diff:	${testInput}/MAprobe.genome.2.1.gaf ${testInput}/MAprobe.genome.3.0.gaf
	cat ${testInput}/MAprobe.genome.2.1.gaf \
	| awk -F'\t' '{ print $$2, $$3, $$4, $$5, $$6, $$7, $$9, $$10, $$12, $$13, $$14, $$15 }' \
        | grep -v ";" |grep -v ":-" |sort > ${scratchDir}/MAprobe.genome.2.1.subset
	cat ${testInput}/MAprobe.genome.3.0.gaf \
	| awk -F'\t' '{ print $$2, $$3, $$4, $$5, $$6, $$7, $$9, $$10, $$12, $$13, $$14, $$15 }' \
        |grep -v ";" |grep -v ":-" |sort > ${scratchDir}/MAprobe.genome.3.0.subset
	diff ${scratchDir}/MAprobe.genome.2.1.subset ${scratchDir}/MAprobe.genome.3.0.subset > $@


${testInput}/MAprobe.genome.3.0.gaf:	${gafDir}/MAprobe.genome.gaf ${testDir}/testMaProbe.txt
	cat ${testDir}/testMaProbe.txt \
	| awk '{ print "grep", $$1, "${gafDir}/MAprobe.genome.gaf"}' \
        | bash > $@

${gafDir}/MAprobe.genome.gaf:	${scratchDir}/MAprobe.genome.raw.gaf ${gafDir}/gene.genome.gaf
	cat $< |sort -k2,2 | scripts/combineFeatures.py > $@

${scratchDir}/MAprobe.genome.raw.gaf:	${inputDir}/MAprobe.genome.bed
	liftOver $< data/GRCh37-lite/hg19.GRCh37-lite.over.chain ${scratchDir}/MAprobe.genome.preGaf.GRCh37-lite.bed /dev/null
	scripts/makeMaProbe.py ${scratchDir}/MAprobe.genome.preGaf.GRCh37-lite.bed $<  > $@

${testInput}/MAprobe.genome.2.1.gaf:	${scratchDir}/MAprobe.genome.gaf21.gaf ${testDir}/testMaProbe.txt 
	cat ${testDir}/testMaProbe.txt \
	| awk '{ print "grep", $$1, "${scratchDir}/MAprobe.genome.gaf21.gaf"}' \
        | bash > $@

#${testDir}/testMaProbe.txt: 
#	hgsql goat -BNe "SELECT FeatureId FROM gaf2_1 WHERE featureType='MAprobe' AND compositeType = 'genome' ORDER BY rand() LIMIT 500" |sort > $@


${scratchDir}/MAprobe.genome.gaf21.gaf:
	zcat ${gaf21File} \
	| awk -F'\t' '$$3 == "MAprobe" && $$9 == "genome" { print }' > $@

#
# When comparing the MAprobes to pre-miRNAs, note that most of the composite
# coordinates WILL be different, because in many cases the pre-miRNA has been
# moved to the other strand by miRBase.  So, omit that from comparison.
#
${testOutput}/MAprobe.pre-miRNA.diff:	${testInput}/MAprobe.pre-miRNA.2.1.gaf ${testInput}/MAprobe.pre-miRNA.3.0.gaf
	cat ${testInput}/MAprobe.pre-miRNA.2.1.gaf \
	| awk -F'\t' '{ print $$2, $$3, $$4, $$5, $$6, $$7, $$8, $$9, $$10, $$13, $$14 }' \
	|sort |sed 's/:-//g' |sed 's/:+//g' \
	> ${scratchDir}/MAprobe.pre-miRNA.2.1.subset
	cat ${testInput}/MAprobe.pre-miRNA.3.0.gaf \
	| awk -F'\t' '{ print $$2, $$3, $$4, $$5, $$6, $$7, $$8, $$9, $$10, $$13, $$14 }' \
	|sort |sed 's/:-//g' |sed 's/:+//g' \
	> ${scratchDir}/MAprobe.pre-miRNA.3.0.subset
	diff ${scratchDir}/MAprobe.pre-miRNA.2.1.subset ${scratchDir}/MAprobe.pre-miRNA.3.0.subset > $@

${testInput}/MAprobe.pre-miRNA.3.0.gaf:     ${gafDir}/MAprobe.pre-miRNA.gaf ${scratchDir}/MAprobe.pre-miRNA.gaf21.gaf 
	- cat ${scratchDir}/MAprobe.pre-miRNA.gaf21.gaf \
	| awk '{ print "grep", $$2, "${gafDir}/MAprobe.pre-miRNA.gaf"}' \
	| bash > ${scratchDir}/MAprobe.pre-miRNA.3.0.superset.gaf
	- cat ${scratchDir}/MAprobe.pre-miRNA.gaf21.gaf \
	| awk -F'\t' '{ print "grep \"" $$8 "\" ${scratchDir}/MAprobe.pre-miRNA.3.0.superset.gaf"}' \
	| bash |sort |uniq > $@	


${testInput}/MAprobe.pre-miRNA.2.1.gaf:     ${testInput}/MAprobe.pre-miRNA.3.0.gaf 
	- cat ${testInput}/MAprobe.pre-miRNA.3.0.gaf \
	| awk '{ print "grep", $$2, "${scratchDir}/MAprobe.pre-miRNA.gaf21.gaf"}' \
	| bash > ${scratchDir}/MAprobe.pre-miRNA.2.1.superset.gaf
	- cat ${testInput}/MAprobe.pre-miRNA.3.0.gaf \
	| awk -F'\t' '{ print "grep \"" $$8 "\" ${scratchDir}/MAprobe.pre-miRNA.2.1.superset.gaf"}' \
	| bash |sort |uniq > $@	


${gafDir}/MAprobe.pre-miRNA.gaf:	${scratchDir}/MAprobe.genome.raw.gaf ${gafDir}/pre-miRNA.genome.gaf ${inputDir}/pre-miRNA.genome.bed
	scripts/gafToBed.py ${gafDir}/pre-miRNA.genome.gaf \
	> ${scratchDir}/pre-miRNA.genome.GRCh37-lite.bed
	liftOver ${scratchDir}/pre-miRNA.genome.GRCh37-lite.bed \
	data/GRCh37-lite/GRCh37-lite.hg19.over.chain \
	${scratchDir}/pre-miRNA.genome.hg19.bed /dev/null
	scripts/maProbeToComposite.py ${scratchDir}/MAprobe.genome.raw.gaf \
	${gafDir}/pre-miRNA.genome.gaf ${scratchDir}/pre-miRNA.genome.hg19.bed > $@



${scratchDir}/MAprobe.pre-miRNA.gaf21.gaf:
	zcat ${gaf21File} \
	| awk -F'\t' '$$3 == "MAprobe" && $$9 == "pre-miRNA" { print }' > $@

#
# When comparing MAprobes to miRNAs, here are the things to NOT look at
# - the miRNA nmes now include miR instead of mir
# - the miRNA names also include suffices (-3p, -5p) where they didn't before.
#   Don't even look at the miRNA name.  The accession is good enough.
# - in any cases, the miRNA coordinates are now on a different strand.  This has two
#   consequences.  First, in cases where the composite (miRNA) does not fully align to
#   the probe, the set of aligned bases will be thrown off.  Second, the composite 
#   coordinates will be in reverse order.  So, don't look at the composite coordinates.
# That doesn't leave much...
#
${testOutput}/MAprobe.miRNA.diff:	${testInput}/MAprobe.miRNA.2.1.gaf ${testInput}/MAprobe.miRNA.3.0.gaf
	cat ${testInput}/MAprobe.miRNA.2.1.gaf \
	| awk -F'\t' '{ split($$8, tokens, "|"); $$8 = tokens[2];  print $$2, $$3, $$4, $$5, $$6, $$7, $$8, $$9, $$10, $$13, $$14, $$16 }' \
	|sort |uniq > ${scratchDir}/MAprobe.miRNA.2.1.subset
	cat ${testInput}/MAprobe.miRNA.3.0.gaf \
	| awk -F'\t' '{ split($$8, tokens, "|"); $$8 = tokens[2]; print $$2, $$3, $$4, $$5, $$6, $$7, $$8, $$9, $$10,  $$13, $$14, $$16 }' \
	|sort |uniq > ${scratchDir}/MAprobe.miRNA.3.0.subset
	diff -i ${scratchDir}/MAprobe.miRNA.2.1.subset ${scratchDir}/MAprobe.miRNA.3.0.subset > $@

${testInput}/MAprobe.miRNA.3.0.gaf:     ${scratchDir}/MAprobe.miRNA.gaf21.gaf ${gafDir}/MAprobe.miRNA.gaf 
	- cat ${scratchDir}/MAprobe.miRNA.gaf21.gaf \
	| awk -F'\t' '{ split($$8, tokens, "|"); print "grep " $$2 ".*" tokens[2] " ${gafDir}/MAprobe.miRNA.gaf"}' \
	| bash |sort |uniq > $@	


${testInput}/MAprobe.miRNA.2.1.gaf:     ${scratchDir}/MAprobe.miRNA.gaf21.gaf ${gafDir}/MAprobe.miRNA.gaf
	- cat ${gafDir}/MAprobe.miRNA.gaf \
	| awk -F'\t' '{ split($$8, tokens, "|"); print "grep " $$2 ".*" tokens[2] " ${scratchDir}/MAprobe.miRNA.gaf21.gaf"}' \
	| bash |sort |uniq > $@	


${gafDir}/MAprobe.miRNA.gaf:	${scratchDir}/MAprobe.genome.raw.gaf ${gafDir}/miRNA.genome.gaf ${inputDir}/miRNA.genome.bed
	scripts/gafToBed.py ${scratchDir}/miRNA.genome.uncombined.gaf \
	> ${scratchDir}/miRNA.genome.GRCh37-lite.bed
	liftOver ${scratchDir}/miRNA.genome.GRCh37-lite.bed \
	data/GRCh37-lite/GRCh37-lite.hg19.over.chain \
	${scratchDir}/miRNA.genome.hg19.bed /dev/null
	scripts/maProbeToComposite.py ${scratchDir}/MAprobe.genome.raw.gaf ${scratchDir}/miRNA.genome.uncombined.gaf  > $@

${scratchDir}/MAprobe.miRNA.gaf21.gaf:
	zcat ${gaf21File} \
	| awk -F'\t' '$$3 == "MAprobe" && $$9 == "miRNA" { print }' > $@
#

${testOutput}/MAprobe.transcript.diff:	${testInput}/MAprobe.transcript.2.1.gaf ${testInput}/MAprobe.transcript.3.0.gaf
	cat ${testInput}/MAprobe.transcript.2.1.gaf \
	| awk -F'\t' '{ print $$2, $$3, $$4, $$5, $$6, $$7, $$8, $$9, $$10, $$11, $$13, $$14, $$15, $$16 }' \
        > ${scratchDir}/MAprobe.transcript.2.1.subset
	cat ${testInput}/MAprobe.transcript.3.0.gaf \
	| awk -F'\t' '{ print $$2, $$3, $$4, $$5, $$6, $$7, $$8, $$9, $$10, $$11, $$13, $$14, $$15, $$16 }' \
        > ${scratchDir}/MAprobe.transcript.3.0.subset
	diff ${scratchDir}/MAprobe.transcript.2.1.subset ${scratchDir}/MAprobe.transcript.3.0.subset > $@

${testInput}/MAprobe.transcript.3.0.gaf:     ${testDir}/testTranscripts.txt ${gafDir}/MAprobe.transcript.gaf 
	cat ${testDir}/testTranscripts.txt \
       | awk '{ print "grep \"" $$1 "\" ${gafDir}/MAprobe.transcript.gaf"}'  \
	| bash |sort -k2,2 > $@       

${testInput}/MAprobe.transcript.2.1.gaf:     ${testDir}/testTranscripts.txt ${scratchDir}/MAprobe.genome.gaf21.gaf                                
	cat ${testDir}/testTranscripts.txt \
       | awk '{ print "grep \"" $$1 "\" ${scratchDir}/MAprobe.transcript.gaf21.gaf"}'\
	| bash |sort -k2,2  > $@       

${gafDir}/MAprobe.transcript.gaf:	${scratchDir}/MAprobe.genome.raw.gaf ${gafDir}/transcript.genome.gaf ${inputDir}/transcript.genome.bed
	scripts/maProbeToComposite.py ${scratchDir}/MAprobe.genome.raw.gaf ${gafDir}/transcript.genome.gaf ${inputDir}/transcript.genome.bed > $@

${scratchDir}/MAprobe.transcript.gaf21.gaf:
	zcat ${gaf21File} \
	| awk -F'\t' '$$3 == "MAprobe" && $$9 == "transcript" { print }' > $@

${scratchDir}/AffySNP.genome.gaf21.gaf: 
	zcat ${gaf21File} \
	| awk -F'\t' '$$3 == "AffySNP" && $$9 == "genome" { print }' > $@

${gafDir}/AffySNP.genome.gaf:	${scratchDir}/AffySNP.genome.gaf21.gaf ${gafDir}/gene.genome.gaf
	cp $< $@

