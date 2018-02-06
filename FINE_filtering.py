from numpy import median

def fine_filtering(templates, reads, k, f):

    number_of_templates = len(templates)
    length_of_read = len(reads[0])

    template_coverage = [0] * number_of_templates


    for ind, template in enumerate(templates):

        extended_template = template + template[:k - 1]
        template_length = len(extended_template)
        all_kmer_templates = [extended_template[i:i + k] for i in range(0, template_length - k + 1)]
        total_number_of_template_kmers = len(all_kmer_templates)

        kmer_templates = set(all_kmer_templates)
        number_of_kmers = len(kmer_templates)

        kmer_matches = [0]*number_of_kmers

        for read in reads:
            kmer_reads = [read[i:i + k] for i in range(0, length_of_read - k + 1)]
            read_kmer_matches = [i for i, item in enumerate(kmer_templates) if item in kmer_reads]

            if len(read_kmer_matches) >= total_number_of_template_kmers - f:
                for match in read_kmer_matches:
                    kmer_matches[match] += 1

        template_coverage[ind] = median(kmer_matches)

    return template_coverage
