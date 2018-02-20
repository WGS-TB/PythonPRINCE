from Bio import SeqIO
from COARSE_filtering import coarse_filtering
from FINE_filtering import fine_filtering




def compute_match_score(listOfGenomes, templates, templateKmers, filteringKmerLength, matchingKmerLength, f0):
    numberOfGenomes = len(listOfGenomes)
    matchScore = [[]] * numberOfGenomes
    genomeCoverage=[]
    for num, genome in enumerate(listOfGenomes):
        try:
            reads1 = list(SeqIO.parse(genome + "1.fq", "fastq"))
            reads2 = list(SeqIO.parse(genome + "2.fq", "fastq"))
        except:
            reads1 = list(SeqIO.parse(genome + "1.fastq", "fastq"))
            reads2 = list(SeqIO.parse(genome + "2.fastq", "fastq"))

        coverage = len(reads1)/1000.0 #assuming all genomes are roughly the same length - # of reads should be proportional to coverage
        genomeCoverage.append(coverage)
        recruitedReads = coarse_filtering(reads1, reads2, filteringKmerLength, templateKmers)
        matchScore[num] = fine_filtering(templates, recruitedReads, matchingKmerLength, f0)

    maxCov = max(genomeCoverage)
    coverage = [x / maxCov for x in genomeCoverage]

    for inx, g in enumerate(matchScore):
        matchScore[inx] = [t / coverage[inx] for t in matchScore[inx]]
    return matchScore


