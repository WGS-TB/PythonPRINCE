from Bio import SeqIO
from COARSE_filtering import coarse_filtering
from FINE_filtering import fine_filtering
from itertools import chain

def check_file_exists(itr8tr):
    first=next(itr8tr)
    return(chain([first],itr8tr))

def compute_match_score(genome, templates, templateKmers, filteringKmerLength, matchingKmerLength, f0):    
    try:
        reads1 = check_file_exists(SeqIO.parse(genome + "1.fq", "fastq"))
        reads2 = check_file_exists(SeqIO.parse(genome + "2.fq", "fastq"))
    except:
        try:
            reads1 = check_file_exists(SeqIO.parse(genome + "1.fastq", "fastq"))
            reads2 = check_file_exists(SeqIO.parse(genome + "2.fastq", "fastq"))
        except:
            try:
                reads1 = check_file_exists(SeqIO.parse(genome + "_1.fq", "fastq"))
                reads2 = check_file_exists(SeqIO.parse(genome + "_2.fq", "fastq"))
            except:
                try:
                    reads1 = check_file_exists(SeqIO.parse(genome + "_1.fastq", "fastq"))
                    reads2 = check_file_exists(SeqIO.parse(genome + "_2.fastq", "fastq"))
                except:
                    raise IOError("Can not open target file %s." % genome)

    #Run reads through Coarse Filtering to drastically reduce computation for Fine Filtering
    nucleotides_seen,recruitedReads = coarse_filtering(chain(reads1,reads2), filteringKmerLength, templateKmers)
    coverage = nucleotides_seen/1000000.0 #assuming all genomes are roughly the same length - nucleotides seen should be proportional to coverage

    #Run reads through Fine Filtering to get score for each template
    matchScore = fine_filtering(templates, recruitedReads, matchingKmerLength,f0)

    #Normalize score by adjusting for coverage
    matchScore = [t/coverage for t in matchScore]    
    return matchScore


