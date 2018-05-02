from Bio import SeqIO
from Kmer_Generator import kmerGenerator
from boost import run_boosts 
from query_sample import test_target
import argparse
import warnings

DEFAULT_K = 9
DEFAULT_BOOST_OUTPUT = "training_data.txt"

def main():
   
    parser = argparse.ArgumentParser(description='Prince Options.')
    
    parser.add_argument('-bo', '--boost_output', default=DEFAULT_BOOST_OUTPUT,
                help="output file for training data / training data used to predict copy numbers for queries")
    parser.add_argument('-to', '--target_output', default="results/predictions.csv",
                help="output file for query copy number predictions")
    parser.add_argument('-tmp','--templates', default="templates.fasta",
                help="VNTR templates. Default is for M.TB")
    parser.add_argument('-tf', '--target_file', default=None,
                help="target genome names in a text file")
    parser.add_argument('-bf', '--boosting_file', default=None,
                help="training genome file names in a text file")
    parser.add_argument('-k', '--k', default=DEFAULT_K,type=int,
                help="Kmer size used during read recruitment.")
    parser.add_argument('-cn', '--copynumber', default=1,type=int,
                help="Copy number for training genome.")

    prince_options = parser.parse_args()

    #Safety check:
    if prince_options.k != DEFAULT_K and prince_options.boost_output == DEFAULT_BOOST_OUTPUT:
	warnings.warn("Warning: Target kmer size does not equal training settings. May lead to inaccurate predictions.")


    #Template data initialized
    templates = list(SeqIO.parse(prince_options.templates, "fasta"))
    templateNames = [t.id for t in templates]
    templates = [str(t.seq) for t in templates]

    #Generate k-mers
    templateKmers = kmerGenerator(templates, prince_options.k)

    if prince_options.boosting_file != None:
        run_boosts(prince_options, templates,templateNames,templateKmers)

    if prince_options.target_file != None:
        test_target(prince_options, templates, templateNames,templateKmers)

if __name__ == '__main__':
    main()
