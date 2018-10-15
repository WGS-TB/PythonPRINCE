from Bio.Seq import Seq
def fine_filtering(templates, reads, k):
    number_of_templates = len(templates)
    template_coverage = [0] * number_of_templates

    for ind, template in enumerate(templates):
        extended_template = template + template[:k - 1]
        reverse_template = Seq(extended_template).reverse_complement()
        template_length = len(extended_template)
        template_kmers = set([extended_template[i:i + k] for i in range(0, template_length - k + 1)] +\
                            [reverse_template[i:i + k] for i in range(0, template_length - k + 1)])

        for read in reads:
            if read in template_kmers:
                template_coverage[ind] += 17

    return template_coverage
