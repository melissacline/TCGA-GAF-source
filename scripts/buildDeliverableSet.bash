#!/usr/bin/env bash
#
# This script is to be run as a final step in splitting the properly-numbered GAF data
# into indivitual feature-component files, and assembling some of those files into
# functional groups.
#
# This script should be run in the GAF deliverables directory.  It takes one argument:
# the completed GAF file.
#

splitSuperset.py $1

cat gene.genome.gaf transcript.genome.gaf transcript.gene.gaf componentExon.genome.gaf \
    componentExon.gene.gaf componentExon.transcript.gaf compositeExon.genome.gaf \
    compositeExon.gene.gaf compositeExon.transcript.gaf junction.genome.gaf \
    junction.transcript.gaf > geneSet.gaf

cat pre-miRNA.genome.gaf miRNA.genome.gaf miRNA.pre-miRNA.gaf > miRnaSet.gaf

cat MAprobe.genome.gaf MAprobe.transcript.gaf MAprobe.pre-miRNA.gaf \
    MAprobe.miRNA.gaf AffySNP.genome.gaf > probeSet.gaf
