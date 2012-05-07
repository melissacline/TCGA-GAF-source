DROP TABLE IF EXISTS knownIsoformsClustersMerged;
CREATE TABLE knownIsoformsClustersMerged (
    clusterId int not null,	# Unique id for transcript cluster (aka gene)
    transcript varchar(255) not null,	# Corresponds to name in knownGene table, transcript in knownCanonical
              #Indices
    INDEX(clusterId),
    UNIQUE(transcript(12))
);
