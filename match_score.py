from Bio import SeqIO
#from modified_COARSE_filtering import coarse_filtering
#from modified_FINE_filtering import fine_filtering
from COARSE_filtering import coarse_filtering
from FINE_filtering import fine_filtering



def compute_match_score(listOfGenomes, templates, templateKmers, filteringKmerLength, matchingKmerLength, f0):
    BOOSTING_LENGTH = 703 #lines in boosting reads (4x number of reads)
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
	print("Done reading")
	reads = reads1+reads2
	reads = [read[n:n+25] for read in reads for n in range(0,len(reads[0])-25+1,25)]
	coverage = len(reads)/1000000.0 #assuming all genomes are roughly the same length - # of reads should be proportional to coverage
        genomeCoverage.append(coverage)
        recruitedReads = coarse_filtering(reads, filteringKmerLength, templateKmers)
        print(len(recruitedReads))
	matchScore[num] = fine_filtering(templates, recruitedReads, matchingKmerLength,f0)
	
        
    for inx, g in enumerate(matchScore):
        matchScore[inx] = [t/genomeCoverage[inx] for t in matchScore[inx]]
    return matchScore


