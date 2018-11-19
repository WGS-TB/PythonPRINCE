from Bio import SeqIO
from Bio.Seq import Seq
from prince.coarse_filtering import coarse_filtering
from prince.fine_filtering import fine_filtering
from prince.kmer_generator import kmer_generator
import itertools
import gzip
import os.path

def combine_records(record1,record2):
    if record2 is None:
        return record1
    else:
        return itertools.chain(record1,record2)

def get_reads_records(filename):
    '''
    Retrieves records for the reads.
    Inputs:
    - (str) filename: a path or prefix to the read file(s)
    Outputs:
    - record1, record2: SeqIO records for the reads
    - gzip_handle1, gzip_handle2: gzip file handles to the reads.
                                  The purpose of these two outputs is so that if the files are contained
                                  in GZIP files, the handles can be closed to make things more efficient
    ''' 
    record1 = None
    record2 = None
    gzip_handle1 = None
    gzip_handle2 = None
    # First, check whether filename is actually just a path
    columns = filename.split('.')
    if columns[-1] == 'gz':
        try:
            gzip_handle1 = gzip.open(filename,'rt')
            record1 = SeqIO.parse(gzip_handle1,'fastq')
        except:
            raise IOError("Cannot open target file %s." % filename)
    elif columns[-1] == 'fastq' or columns[-1] == 'fq':
        try:
                record1 = SeqIO.parse(filename,'fastq')
        except:
            raise IOError("Can not open target file %s." % filename)
    # Otherwise, check all possible file extensions
    else:
        delimiters = ['','_']
        extensions = ['','.fastq','.fq']
        gzipped = ['','.gz']
        # Check whether the reads is contained only in a single file:
        for delim,exten,gz in itertools.product(delimiters,extensions,gzipped):
            path = filename + delim + exten + gz
            if os.path.exists(path):
                if gz == '.gz':
                    try:
                        gzip_handle1 = gzip.open(path,'rt')
                        record1 = SeqIO.parse(gzip_handle1,'fastq')
                    except:
                        raise IOError("Cannot open target file " + filename)
                else:
                    try:
                        record1 = SeqIO.parse(path,'fastq')
                    except:
                        raise IOError("Cannot open target file " + filename)
        # Otherwise, check whether the reads are paired end
        if record1 is None:
            for delim,exten,gz in itertools.product(delimiters,extensions,gzipped):
                path1 = filename + delim + "1" + exten + gz
                path2 = filename + delim + "2" + exten + gz
                if os.path.exists(path1) and os.path.exists(path2):
                    if gz == '.gz':
                        try:
                            gzip_handle1 = gzip.open(path1,'rt') 
                            gzip_handle2 = gzip.open(path2,'rt')
                            record1 = SeqIO.parse(gzip_handle1,'fastq') 
                            record2 = SeqIO.parse(gzip_handle2,'fastq')
                        except:
                            raise IOError("Cannot open target file " + filename)
                    else:
                        try:
                            record1 = SeqIO.parse(path1,'fastq')
                            record2 = SeqIO.parse(path2,'fastq')
                        except:
                            raise IOError("Cannot open target file " + filename)
        if record1 is None:
            raise IOError("Cannot open target file " + filename)
    return record1,record2,gzip_handle1,gzip_handle2
                        
            
def compute_match_score(filename, template_obj, kmerLength, primers):
    '''
    Inputs:
    - (str) data_prefix: the prefix of the NGS dataset paths
    '''
    record1, record2, gzip1, gzip2 = get_reads_records(filename) 
    #Run reads through Coarse Filtering to drastically reduce computation for Fine Filtering
    reads = combine_records(record1,record2)
    
    
    flanking_sequences = set(kmer for kmers in kmer_generator([sequence for primer_set in primers.values() for sequence in primer_set], kmerLength, extension=False).values() for kmer in kmers)
    
    nucleotides_seen,recruitedReads = coarse_filtering(reads, kmerLength, template_obj["Kmers"], flanking_sequences)
    
    # Close all the things
    record1.close()
    if record2 is not None:
        record2.close()
    if gzip1 is not None:
        gzip1.close()
    if gzip2 is not None:
        gzip2.close()

    coverage = nucleotides_seen/1000000.0 #assuming all genomes are roughly the same length - nucleotides seen should be proportional to coverage

    #Run reads through Fine Filtering to get score for each template
    matchScore, flanking_coverage = fine_filtering(template_obj, recruitedReads, kmerLength, primers)
    print(matchScore)
    print(flanking_coverage)
    #Normalize score by adjusting for coverage and flanking coverage
    #matchScore = [t/coverage for t in matchScore]
    #matchScore = [matchScore[i]/(coverage*(1+flanking_coverage[i])) for i in range(len(matchScore))]
    matchScore = [score/float(1+flanking_coverage[i]) for score in matchScore]
    print(matchScore)
    print("\n")
    return matchScore

