from Bio import SeqIO
from Kmer_Generator import kmerGenerator
from boost import run_boosts 
from query_sample import test_target
import argparse

import time
start_time = time.time()

#To Train: need: -bg -b -bo -cn --- optional: -c -m -f
#To Query: need: -to -t  	--- optional: -c -m -f -bo


def main():
   
    parser = argparse.ArgumentParser(description='get matchscores.')
    
    parser.add_argument('-bo', '--boost_output', default="output/matchscores.csv",
                help="output file for matchscores")
    parser.add_argument('-to', '--target_output', default="results/predictions.csv",
                help="output file for query copy number predictions")
    parser.add_argument('-tmp','--templates', default="templates.fasta",
                help="VNTR templates")
    parser.add_argument('-tf', '--target_file', default=None,
                help="target genome names in a text file")
    parser.add_argument('-bf', '--boosting_file', default=None, #GET RID OF LATER
                help="boosting genome file names in a text file")
    parser.add_argument('-c', '--coarse', default=11,type=int,
                help="coarse filtering kmer length")
    parser.add_argument('-f', '--fine', default=9,type=int,
                help="fine filtering kmer length")
    parser.add_argument('-s', '--screen', default=10,type=int,
                help="Kmer mismatches allowed during fine filtering.")
    parser.add_argument('-cn', '--copynumber', default=1,type=int,
                help="Copy number for training genome.")

    prince_options = parser.parse_args()


    # Variable declarations and imports
    # # #####################################

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
