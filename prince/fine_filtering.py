def fine_filtering(templates, reads, k):
    number_of_templates = len(templates)
    length_of_read = len(reads[0])
    template_coverage = [0] * number_of_templates

    for ind, template in enumerate(templates):
        extended_template = template + template[:k - 1]
        template_length = len(extended_template)
        template_kmers = set([extended_template[i:i + k] for i in range(0, template_length - k + 1)])

        for read in reads:
            read_kmers = (read[i:i + k] for i in range(0, length_of_read - k + 1))
            num_read_matches = sum([1 for kmer in read_kmers if kmer in template_kmers])

	    if num_read_matches >= length_of_read-k+1:
                template_coverage[ind] += num_read_matches
    return template_coverage