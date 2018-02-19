from Bio import SeqIO
import numpy as np
from sklearn.cluster import KMeans
from Kmer_Generator import kmerGenerator
from COARSE_filtering import coarse_filtering
from FINE_filtering import fine_filtering
from match_score import compute_match_score
import argparse

parser = argparse.ArgumentParser(description='get matchscores.')
parser.add_argument('-tg', '--target', default=None,
            help="target genome reads")
parser.add_argument('-c', '--coarse', default=25,
            help="coarse filtering kmer length")
parser.add_argument('-m', '--match', default=15,
            help="fine filtering kmer length")
parser.add_argument('-f', '--filter', default=10,
            help="Kmer mismatches allowed during fine filtering.")
opts = parser.parse_args()


# Variable declarations and imports
# # #####################################
filteringKmerLength = opts.coarse
matchingKmerLength = opts.match


f0 = opts.filter

numberOfClusters = 9

templates = list(SeqIO.parse("templates.fasta", "fasta"))
templates = [str(t.seq) for t in templates]



#Insert names of TRAINING genomes here. No extentions, just the name.
boostingGenomeNames = ["kurono-sequence1_150bp_15x"]



#Insert names of YOUR genomes here. No extentions, just the name.
targetGenomeNames = [opts.target]


# Determining the match score
# # #####################################


#Generate k-mers
templateKmers = kmerGenerator(templates, filteringKmerLength)




#FINE FILTERING, COARSE FILTERING + MATCH SCORE for TRAINING genomes
boostingMatchScore = compute_match_score(boostingGenomeNames, templates, templateKmers,
                                         filteringKmerLength, matchingKmerLength, f0)



#FINE FILTERING, COARSE FILTERING + MATCH SCORE for YOUR genomes
targetMatchScore = compute_match_score(targetGenomeNames, templates, templateKmers,
                                       filteringKmerLength, matchingKmerLength, f0)

######################
#Clustering

for index, template in enumerate(templates):
    # X = [9, 11, 21, 22, 27, 29, 31, 40, 53, 49, 62, 69]
    X = [row[index] for row in boostingMatchScore]
    # X = [j for i in boostingMatchScore for j in i]
    print X
    predict = [row[index] for row in targetMatchScore]
    # predict = [j for i in targetMatchScore for j in i]

    X = (np.array(X)).reshape(-1, 1)
    predict = (np.array(predict)).reshape(-1, 1)


    kmeans = KMeans(n_clusters=numberOfClusters, random_state=0).fit(X)

    output = kmeans.predict(predict)
    print kmeans.cluster_centers_




    # # https://stackoverflow.com/questions/44888415/how-to-set-k-means-clustering-labels-from-highest-to-lowest-with-python
    # The values in the lookuptable correspond to the values of the CNVs (The order does)

    idx = np.argsort(kmeans.cluster_centers_.sum(axis=1))
    lookupTable = np.zeros_like(idx)
    lookupTable[idx] = np.arange(numberOfClusters)
    # print idx
    # print "lookup table", lookupTable
    print lookupTable[output.labels_]
    #



    #run k-means on it
    #sort centers of clusters
     #check what sorted center of cluster target belongs to



