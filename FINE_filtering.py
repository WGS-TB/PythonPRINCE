from numpy import median

def fine_filtering(templates, reads, k, f=15):
    number_of_templates = len(templates)
    length_of_read = len(reads[0])
#    reads = [read[n:n+25] for read in reads for n in range(length_of_read-25+1)] #to account for different read lengths
    template_coverage = [0] * number_of_templates
    '''
    for ind, template in enumerate(templates):
        extended_template = template + template[:k - 1]
        template_length = len(extended_template)
        kmer_templates = [extended_template[i:i + k] for i in range(0, template_length - k + 1)]
        total_number_of_template_kmers = len(kmer_templates)

        number_of_kmers = len(kmer_templates)
        kmer_matches = [0]*number_of_kmers

        for read in reads:
            kmer_reads = set([read[i:i + k] for i in range(0, length_of_read - k + 1)])

            read_kmer_matches = [i for i, item in enumerate(kmer_templates) if item in kmer_reads]
            if len(read_kmer_matches) >= total_number_of_template_kmers-f:
                for match in read_kmer_matches:
                    kmer_matches[match] += 1
                    # print kmer_matches

        template_coverage[ind] = median(kmer_matches)
    '''
    for ind, template in enumerate(templates):
        extended_template = template + template[:k - 1]
        template_length = len(extended_template)
        template_kmers = set([extended_template[i:i + k] for i in range(0, template_length - k + 1)])
        length_of_read = len(reads[0])

        for read in reads:
            read_kmers = [read[i:i + k] for i in range(0, length_of_read - k + 1)]


            read_matches = [1 for i, item in enumerate(read_kmers) if item in template_kmers]
            if len(read_matches) >= length_of_read-f:
                template_coverage[ind] += len(read_matches)
                    # print kmer_matches
    return template_coverage


