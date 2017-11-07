from Bio import SeqIO
from COARSE_filtering import coarse_filtering
from FINE_filtering import fine_filtering




def compute_match_score(listOfGenomes,genomeCoverage, templates, templateKmers, filteringKmerLength, matchingKmerLength, f0):
    numberOfGenomes = len(listOfGenomes)
    matchScore = [[]] * numberOfGenomes

    for num, genome in enumerate(listOfGenomes):

        reads1 = list(SeqIO.parse(genome + "_1.fq", "fastq"))
        reads2 = list(SeqIO.parse(genome + "_2.fq", "fastq"))

        recruitedReads = coarse_filtering(reads1, reads2, filteringKmerLength, templateKmers)
        matchScore[num] = fine_filtering(templates, recruitedReads, matchingKmerLength, f0)

    maxCov = float(max(genomeCoverage))
    coverage = [x / maxCov for x in genomeCoverage]

    for inx, g in enumerate(matchScore):
        matchScore[inx] = [t / coverage[inx] for t in matchScore[inx]]
    return matchScore


