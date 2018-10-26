from Bio import SeqIO
from prince.coarse_filtering import coarse_filtering
from prince.fine_filtering import fine_filtering
from itertools import chain
import gzip

def check_file_exists(itr8tr):
    first=next(itr8tr)
    return(chain([first],itr8tr))

def compute_match_score(genome, genome_reverse, templates, templateKmers, kmerLength):
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
                        with gzip.open(genome + "_1.fastq.gz", "rt") as handle:
                            reads1 = check_file_exists(SeqIO.parse(handle, "fastq"))

                        with gzip.open(genome + "_2.fastq.gz", "rt") as handle:
                            reads2 = check_file_exists(SeqIO.parse(handle, "fastq"))
                    except:
                        try:
                            with gzip.open(genome + "_1.fq.gz", "rt") as handle:
                                reads1 = check_file_exists(SeqIO.parse(handle, "fastq"))

                            with gzip.open(genome + "_2.fq.gz", "rt") as handle:
                                reads2 = check_file_exists(SeqIO.parse(handle, "fastq"))
                        except:
                            try:
                                if genome_reverse == "":
                                    reads1 = check_file_exists(SeqIO.parse(genome, "fastq"))
                                    reads2 = iter(())
                                else:
                                    reads1 = check_file_exists(SeqIO.parse(genome, "fastq"))
                                    reads2 = check_file_exists(SeqIO.parse(genome_reverse, "fastq"))
                            except:
                                try:
                                    if genome_reverse == "":
                                        with gzip.open(genome, "rt") as handle:
                                            reads1 = check_file_exists(SeqIO.parse(handle, "fastq"))
                                        reads2 = iter(())
                                    else:
                                        with gzip.open(genome, "rt") as handle:
                                            reads1 = check_file_exists(SeqIO.parse(handle, "fastq"))
                                        with gzip.open(genome_reverse, "rt") as handle:
                                            reads2 = check_file_exists(SeqIO.parse(handle, "fastq"))
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

