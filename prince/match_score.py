from Bio import SeqIO
from prince.coarse_filtering import coarse_filtering
from prince.fine_filtering import fine_filtering
from itertools import chain
import gzip
import os.path

def check_file_exists(itr8tr):
    first=next(itr8tr)
    return(chain([first],itr8tr))

def compute_match_score(genome, templates, templateKmers, kmerLength):
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
                    try:
                        handle1 = gzip.open(genome + "_1.fastq.gz", "rt")
                        reads1 = check_file_exists(SeqIO.parse(handle1, "fastq"))

                        handle2 = gzip.open(genome + "_2.fastq.gz", "rt")
                        reads2 = check_file_exists(SeqIO.parse(handle2, "fastq"))
                    except:
                        try:
                            reads1 = check_file_exists(SeqIO.parse(genome, "fastq"))
                            reads2 = iter(())
                        except:
                            raise IOError("Can not open target file %s." % genome)

    #Run reads through Coarse Filtering to drastically reduce computation for Fine Filtering
    nucleotides_seen,recruitedReads = coarse_filtering(chain(reads1,reads2), kmerLength, templateKmers)
    coverage = nucleotides_seen/1000000.0 #assuming all genomes are roughly the same length - nucleotides seen should be proportional to coverage

    #Run reads through Fine Filtering to get score for each template
    matchScore = fine_filtering(templates, recruitedReads, kmerLength)

    #Normalize score by adjusting for coverage
    matchScore = [t/coverage for t in matchScore]
    return matchScore

