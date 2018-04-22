from Bio import SeqIO
from COARSE_filtering import coarse_filtering
from FINE_filtering import fine_filtering
from itertools import chain


def compute_match_score(genome, templates, templateKmers, filteringKmerLength, matchingKmerLength, f0):
    import time
    start_time = time.time()
    
    try:
        reads1 = SeqIO.parse(genome + "1.fq", "fastq")
        reads2 = SeqIO.parse(genome + "2.fq", "fastq")
    except:
        reads1 = SeqIO.parse(genome + "1.fastq", "fastq")
        reads2 = SeqIO.parse(genome + "2.fastq", "fastq")

    length_of_reads,recruitedReads = coarse_filtering(chain(reads1,reads2), filteringKmerLength, templateKmers)
    print(length_of_reads)
    coverage = length_of_reads/1000000.0 #assuming all genomes are roughly the same length - # of reads should be proportional to coverage

    print(len(recruitedReads))
    print(time.time()-start_time)
    matchScore = fine_filtering(templates, recruitedReads, matchingKmerLength,f0)

    matchScore = [t/coverage for t in matchScore]
    
    return matchScore


