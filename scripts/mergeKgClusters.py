#!/usr/bin/env python

import argparse
import MySQLdb
import MySQLdb.cursors
import re

def transcriptsOverlap(transcriptId1, transcriptId2, cursor):
    cursor.execute("select * from knownGene where name = '%s'" % (transcriptId1))
    assert(cursor.rowcount == 1)
    trans1 = cursor.fetchone()
    cursor.execute("select * from knownGene where name = '%s'" % (transcriptId2))
    assert(cursor.rowcount == 1)
    trans2 = cursor.fetchone()
    if trans1['chrom'] == trans2['chrom'] \
           and trans1['strand'] == trans2['strand']:
        if trans1['exonCount'] == 1 and trans2['exonCount'] == 1:
            #
            # Neither transcript is spliced.  If their transcripts overlap, t
            # then we're done.
            if trans1['txStart'] < trans2['txEnd'] \
                   and trans2['txStart'] < trans1['txEnd']:
                print trans1['name'], "and", trans2['name'], "are unspliced and overlap"
                return(True)
        else:
            #
            # Pull out the exon endpoints, removing trailing commas as needed.
            # If the transcripts share some start or end splice site, we're done.
            exonStarts1 = re.sub(",$", "", trans1['exonStarts']).split(',')
            exonStarts2 = re.sub(",$", "", trans2['exonStarts']).split(',')
            exonEnds1 = re.sub(",$", "", trans1['exonEnds']).split(',')
            exonEnds2 = re.sub(",$", "", trans2['exonEnds']).split(',')
            #
            # First test the 3' splice sites
            for ii in range(1, len(exonStarts1)):
                for jj in range(1, len(exonStarts2)):
                    if exonStarts1[ii] == exonStarts2[jj]:
                        print trans1['name'], "and", trans2['name'], \
                              "share an exonStart", exonStarts1[ii], exonStarts2[jj]
                        return(True)
            #
            # Next test the 5' splice sites
            for ii in range(0, len(exonEnds1) - 1):
                for jj in range(0, len(exonEnds2) - 1):
                    if exonEnds1[ii] == exonEnds2[jj]:
                        print trans1['name'], "and", trans2['name'], \
                              "share an exonEnd", exonEnds1[ii], exonEnds2[jj]
                        return(True)
    return(False)


def clustersOverlap(clusterId1, clusterId2, cursor):
    """Determine if two clusters overlap and should be merged """
    #
    # Step 1: determine if the clusters share any transcript(s) with the
    # same gene symbol
    cursor.execute("""select ki1.transcript as transcript1,
                             ki2.transcript as transcript2
                        from knownIsoforms ki1, kgXref kx1,
                             knownIsoforms ki2, kgXref kx2
                        where ki1.transcript = kx1.kgId and ki1.clusterId = %s
                          and ki2.transcript = kx2.kgId and ki2.clusterId = %s
                          and kx1.geneSymbol = kx2.geneSymbol""" \
                   % (clusterId1, clusterId2))
    for row in cursor.fetchall():
        #
        # Step 2: determine if these transcripts have some overlapping block,
        # and share a splice site if either is spliced.  Also verify that they're
        # on the same strand.
        if transcriptsOverlap(row["transcript1"], row["transcript2"], cursor):
            print "found overlap between", row["transcript1"], "of", clusterId1, \
                  "and", row["transcript2"], "of", clusterId2
            return(True)
    return(False)


def findOverlappingClusters(clusterSet, clusterId1, cursor):
    """Find any clusters that should be merged with the indicated cluster,
    given that they have some transcript with the same symbol and an
    overlapping exon, and they share a splice site if they are spliced.
    Add the indicated clusterId to the cluster set.  If any overlapping clusters
    are found, add them to the cluster set, and then recurse to find any
    clusters that should be merged with them"""
    clusterSet.add(clusterId1)
    cursor.execute("""select clusterId2 from kgOverlappingClusters
                       where clusterId1 = %s""" % (clusterId1))
    for row in cursor.fetchall():
        clusterId2 = row["clusterId2"]
        print "Exploring candidate match betweeen", clusterId1, "and", clusterId2
        if clustersOverlap(clusterId1, clusterId2, cursor):
            print clusterId1, "and", clusterId2, "overlap, clusterSet now", clusterSet, "recursing with", clusterId2
            clusterSet = findOverlappingClusters(clusterSet, clusterId2, cursor)
    return(clusterSet)


def addTranscriptsToCluster(clusterToAdd, clusterToAddTo, table, cursor):
    """Add the transcripts from the first cluster into the second cluster"""
    cursor.execute("""select transcript from knownIsoforms 
                       where clusterId = %s""" % (clusterToAdd))
    for row in cursor.fetchall():
        cursor.execute("""insert into %s
                          (clusterId, transcript) values (%s, '%s')""" \
                       % (table, clusterToAddTo, row["transcript"]))
        
parser = argparse.ArgumentParser()
parser.add_argument('db', type=str, help="database with knownGene table",
                    default="hg19")
parser.add_argument('destinationTable', type=str, default="hg19",
                    help="Table that will contain the merged clusters")
args = parser.parse_args()

db = MySQLdb.connect(read_default_file="~/.my.hg19.cnf")
cursor = db.cursor(MySQLdb.cursors.DictCursor)

cursor.execute("""drop table if exists kgOverlappingClusters""")
cursor.execute("""create table kgOverlappingClusters as
                 select kc1.clusterId as clusterId1, kc2.clusterId as clusterId2
                    from knownCanonical kc1,knownGene kg1, knownCanonical kc2,
                         knownGene kg2
                   where kg1.name = kc1.transcript and kg2.name = kc2.transcript
                     and kc1.clusterId != kc2.clusterId and kc1.chrom = kc2.chrom
                     and kg1.strand = kg2.strand and kg1.txStart < kg2.txEnd
                     and kg2.txStart < kg1.txEnd
                     and kc1.clusterId < kc2.clusterId""")
clustersCopied = set()
cursor.execute("select clusterId from knownCanonical order by clusterId")
for row in cursor.fetchall():
    clusterId = row["clusterId"]
    if clusterId not in clustersCopied:
        print "looking closer at", clusterId
        overlappingClusters = findOverlappingClusters(set(), clusterId, cursor)
        for clusterInOverlap in overlappingClusters:
            print "copying transcripts from", clusterInOverlap, "to", clusterId
            addTranscriptsToCluster(clusterInOverlap, clusterId, 
                                    args.destinationTable, cursor)
            clustersCopied.add(clusterInOverlap)

