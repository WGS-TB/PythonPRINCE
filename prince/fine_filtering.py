from Bio.Seq import Seq
def fine_filtering(template_obj, reads, k, primers):
    templates = template_obj["Sequences"]
    number_of_templates = len(templates)
    template_coverage = [0] * number_of_templates
    primer_coverage = [0] * number_of_templates
    template_names = template_obj["Names"]
    
    for ind, template in enumerate(templates):
        template_name = template_names[ind]
        
        primer_sequences = set([kmer for primer in primers[template_name]\
                            for direction in (primer,str(Seq(primer).reverse_complement()))\
                            for kmer in [direction[i:i + k] for i in range(0, len(direction) - k + 1)]])
        
        template_kmers = set(template_obj["Kmers"][ind])
        
        for read in reads:
            if read in template_kmers:
                template_coverage[ind] += 1
            if read in primer_sequences:
                primer_coverage[ind] += 1
        
    return (template_coverage, primer_coverage)
