from numpy import median

def fine_filtering(templates, reads, k, f):

    number_of_templates = len(templates)
    length_of_read = len(reads[0])
    template_coverage = [0] * number_of_templates

    for ind, template in enumerate(templates):

        template_length = len(template)

        kmer_templates = list(set([template[i:i + k] for i in range(0, template_length - k + 1)]))

        number_of_kmers = len(kmer_templates)

        kmer_matches = [0]*number_of_kmers

        for read in reads:

            kmer_reads = [read[i:i + k] for i in range(0, length_of_read - k + 1)]

            read_kmer_matches = [kmer_templates.index(i) for i in kmer_reads if i in kmer_templates]

            if len(read_kmer_matches) >= f:
                for match in read_kmer_matches:
                    kmer_matches[match] += 1

        template_coverage[ind] = median(kmer_matches)

    return template_coverage


#
#
#
#         # print kmer_templates
#
# print( fine_filtering(["12345867"], ["12345", "23999", "56799"], 2, 0))
#
