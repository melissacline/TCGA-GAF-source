inputDir = data/gaf4_0-init
gafDir = data/gaf4_0
outputBedDir = data/gaf4_0-bed
scratchDir = data/tmp
version = v4_0


# This stops make from deleting intermediate files
.SECONDARY:

# Genes, transcripts, exons, junctions

geneSetGaf = ${gafDir}/gene.genome.${version}.gaf ${gafDir}/transcript.genome.${version}.gaf \
        ${gafDir}/transcript.gene.${version}.gaf \
        ${gafDir}/exon.genome.${version}.gaf ${gafDir}/exon.gene.${version}.gaf \
        ${gafDir}/exon.transcript.${version}.gaf \
        ${gafDir}/junction.genome.${version}.gaf ${gafDir}/junction.gene.${version}.gaf \
        ${gafDir}/junction.transcript.${version}.gaf

# SNP and probesets
snpSetGaf = ${gafDir}/dbSNP.genome.${version}.gaf

MAprobeGeneGaf = ${gafDir}/MAprobe.genome.${version}.gaf
MAprobeTranscriptGaf = ${gafDir}/MAprobe.transcript.${version}.gaf

affySnpGaf = ${gafDir}/AffySNP.genome.${version}.gaf

probeSetGaf = ${MAprobeGeneGaf} ${MAprobeTranscriptGaf} ${affySnpGaf}

# supersets
superSets = ${gafDir}/geneSet.${version}.gaf ${gafDir}/probeSet.$version.gaf

# other inputs
refseqData = ${inputDir}/gencode.v19.metadata.RefSeq
gencodeSet = ${inputDir}/gencode.v19.annotation.gtf
GRCh37.hg19.chain = ${inputDir}/GRCh37-lite.hg19.over.chain
hg19.GRCh37.chain = ${inputDir}/hg19.GRCh37-lite.over.chain
dbSNP_2.1 = ${inputDir}/GAF2.1/TCGA.hg19.June2011.with_dbSNP.gaf.gz


# bigbed files
geneSetBB = ${outputBedDir}/gene.genome.$version.bb ${outputBedDir}/transcript.genome.$version.bb \
        ${outputBedDir}/compositeExon.genome.$version.bb ${outputBedDir}/componentExon.genome.$version.bb \
        ${outputBedDir}/junction.genome.$version.bb
snpSetBB = ${outputBedDir}/dbSNP.genome.$version.bb ${outputBedDir}/AffySNP.genome.$version.bb
probeSetBB = ${outputBedDir}/MAprobe.genome.$version.bb

geneSet = ${geneSetGaf} ${snpSetGaf} ${MAprobeSetGaf} ${affySnpGaf}
bedSet = ${geneSetBB} ${snpSetBB} ${probeSetBB}

all:    ${geneSet} ${superSet} ${bedSet} ${probeSetGaf}

# gtfToGAF.py creates all geneSetGaf files
${geneSetGaf}: ${scratchDir}/gaf.done
	

${scratchDir}/gaf.done: ${gencodeSet} ${refseqData}
	mkdir -p $(dir $@)
	mkdir data/gaf4_0
	scripts/gtfToGAF.py ${gencodeSet} ${refseqData} ${gafDir} ${version}
	touch $@

# The other scripts all query the gafGeneXref table in hg19, so create it
${scratchDir}/gafGeneXref.table:	${gafDir}/gene.genome.${version}.gaf
	mkdir -p $(dir $@)
	scripts/geneToGeneXref_fake.py $< > $@
	echo "delete from gafGeneXref" | hgsql hg19
	hgsqlimport --local hg19 $@


# AffySNPs
# makeAffySnp.py uses the version 2.1 gaf file to get GAF details

${scratchDir}/AffySNP.genome.gaf21.gaf:		${dbSNP_2.1} ${scratchDir}/gafGeneXref.table
	mkdir -p $(dir $@)
	zcat ${dbSNP_2.1} | \
	awk -F'\t' '$$3 == "AffySNP" && $$9 == "genome" { print }' > $@

${affySnpGaf}:  ${inputDir}/AffySNP.genome.bed ${hg19.GRCh37.chain} ${scratchDir}/AffySNP.genome.gaf21.gaf
	mkdir -p $(dir $@)
	liftOver ${inputDir}/AffySNP.genome.bed ${hg19.GRCh37.chain} \
	${scratchDir}/AffySNP.genome.preGaf.GRCh37-lite.bed /dev/null
	count=$$(wc -l ${gafDir}/*.gaf | grep total | awk '{print $$1}'); let count=$$count+1; \
	scripts/makeAffySnp.py -n $$count ${scratchDir}/AffySNP.genome.preGaf.GRCh37-lite.bed \
	${scratchDir}/AffySNP.genome.gaf21.gaf > $@

# MAprobes
${scratchDir}/MAprobe.genome.bed:
	mkdir -p $(dir $@)
	hgsql hg19 -Ne "select * from maProbe" | cut -f2- > $@

${MAprobeGeneGaf}:	${scratchDir}/MAprobe.genome.bed ${hg19.GRCh37.chain} ${scratchDir}/gafGeneXref.table
	mkdir -p $(dir $@)
	liftOver ${scratchDir}/MAprobe.genome.bed ${hg19.GRCh37.chain} \
	${scratchDir}/MAprobe.genome.preGaf.GRCh37-lite.bed /dev/null
	# makeMaProbe.py queries the gafGeneXref table
	scripts/makeMaProbe.py ${scratchDir}/MAprobe.genome.preGaf.GRCh37-lite.bed  > ${scratchDir}/MAprobe.genome.raw.gaf
	count=$$(wc -l ${gafDir}/*.gaf | grep total | awk '{print $$1}'); let count=$$count+1; \
	cat ${scratchDir}/MAprobe.genome.raw.gaf |sort -k2,2 | scripts/combineFeatures.py |\
	awk -v i=$$count 'BEGIN {OFS="\t";FS="\t"}; /.*/{$$1=""; printf "%d\t% s\n",i,$$0; i++}' |\
	cut -f1,3- > $@

${MAprobeTranscriptGaf}:	${MAprobeGeneGaf} ${gafDir}/transcript.genome.${version}.gaf ${scratchDir}/gafGeneXref.table
	mkdir -p $(dir $@)
	count=$$(wc -l ${gafDir}/*.gaf | grep total | awk '{print $$1}'); let count=$$count+1; \
	scripts/maProbeToComposite.py -n $$count ${scratchDir}/MAprobe.genome.raw.gaf ${gafDir}/transcript.genome.${version}.gaf > $@


# dbSNP. This takes several days, use with caution.
${scratchDir}/snp138.hg19.bed:
	mkdir -p $(dir $@)
	hgsql hg19 -Ne "select * from snp138"  | cut -f2-7 | \
	awk 'BEGIN{OFS="\t";}{print $$1,$$2,$$3,$$4,$$5,$$6,$$2,$$3,"0","1",$$3-$$2,0}' > $@


### makeDbSnp.py takes a few days so use with caution
### Using awk to put the unique entryNumber in the first field
### The first number is the linecount of all gaf files created so far, plus one.
${snpSetGaf}:      ${scratchDir}/dbSNP.raw.gaf
	count=$$(wc -l ${gafDir}/*.gaf | grep total | awk '{print $$1}'); let count=$$count+1; \
	cat $< | sort -k2,2 | scripts/combineSnps.py | \
	awk -v i=$$count 'BEGIN{OFS="\t";FS="\t"}; /.*/{$$1=""; printf "%d\t% s\n",i,$$0; i++}' |\
	cut -f1,3- > $@

${scratchDir}/dbSNP.raw.gaf:      ${scratchDir}/snp138.hg19.bed ${hg19.GRCh37.chain} ${scratchDir}/gafGeneXref.table
	liftOver ${scratchDir}/snp138.hg19.bed ${hg19.GRCh37.chain} ${scratchDir}/dbSNP.genome.preGaf.GRCh37-lite.bed /dev/null
	scripts/makeDbSnp.py ${scratchDir}/dbSNP.genome.preGaf.GRCh37-lite.bed ${scratchDir}/snp138.hg19.bed > $@

	count=$$(wc -l ${gafDir}/*.gaf | grep total | awk '{print $$1}'); let count=$$count+1; \
	cat ${scratchDir}/dbSNP.raw.gaf | sort -k2,2 | scripts/combineSnps.py | \
	awk -v i=$$count 'BEGIN{OFS="\t";FS="\t"}; /.*/{$$1=""; printf "%d\t% s\n",i,$$0; i++}' |\
	cut -f1,3- > $@


# create superset gafs
${gafDir}/geneSet.${version}.gaf:	${geneSetGaf}
	cat $< > $@
${gafDir}/probeSet.${version}.gaf:	${probeSetGaf}
	cat $< > $@


# Bed files are created using a shell script that deals with creating thickStart/thickEnd for transcripts
${bedSet}:	 ${scratchDir}/bed.done

${scratchDir}/bed.done:    ${GRCh37.hg19.chain} ${geneSetGaf} ${probeSetGaf} ${snpSetGaf}
	mkdir -p $(dir $@)
	mkdir -p data/gaf4_0-bed
	hgsql hg19 -Ne "select chrom, size from chromInfo" > ${scratchDir}/chrom.sizes
	scripts/gafToBigBed.sh ${gafDir} ${version} ${scratchDir} ${GRCh37.hg19.chain} ${outputBedDir} ${scratchDir}/chrom.sizes
	touch $@


cleanTrash:
	rm -f ${scratchDir}/*
	rmdir ${scratchDir} 2>/dev/null
clean:
	rm -f ${gafDir}/* ${scratchDir}/* ${outputBedDir}/*
	rmdir ${gafDir} ${scratchDir} ${outputBedDir} 2>/dev/null
