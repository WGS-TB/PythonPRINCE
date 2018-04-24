from Bio import SeqIO
from Kmer_Generator import kmerGenerator
from boost import run_boosts 
from query_sample import test_target
import argparse
import warnings

#To Train: need: -bf -bo -cn	--- optional: -c -m -f
#To Query: need: -tf -to  	--- optional: -c -m -f -bo

DEFAULT_COARSE = 9
DEFAULT_FINE = 9
DEFAULT_SCREEN = 8
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
    parser.add_argument('-c', '--coarse', default=DEFAULT_COARSE,type=int,
                help="coarse filtering kmer length")
    parser.add_argument('-f', '--fine', default=DEFAULT_FINE,type=int,
                help="fine filtering kmer length")
    parser.add_argument('-s', '--screen', default=DEFAULT_SCREEN,type=int,
                help="Kmer mismatches allowed during fine filtering.")
    parser.add_argument('-cn', '--copynumber', default=1,type=int,
                help="Copy number for training genome.")

    prince_options = parser.parse_args()

    #Safety check:
    if (prince_options.coarse != DEFAULT_COARSE or prince_options.fine != DEFAULT_FINE or prince_options.screen != DEFAULT_SCREEN)\
	and prince_options.boost_output == DEFAULT_BOOST_OUTPUT:
	warnings.warn("Warning: Target settings (-c or -f or -s) do not equal training settings. May lead to inaccurate predictions.")


    #Template data initialized
    templates = list(SeqIO.parse(prince_options.templates, "fasta"))
    templateNames = [t.id for t in templates]
    templates = [str(t.seq) for t in templates]

    #Generate k-mers
    templateKmers = kmerGenerator(templates, prince_options.coarse)

    if prince_options.boosting_file != None:
        run_boosts(prince_options, templates,templateNames,templateKmers)

    if prince_options.target_file != None:
        test_target(prince_options, templates, templateKmers)

if __name__ == '__main__':
    main()
