from Bio import SeqIO
from Kmer_Generator import kmerGenerator
from boost import run_boosts 
from query_sample import test_target
import argparse


parser = argparse.ArgumentParser(description='get matchscores.')
parser.add_argument('-t', '--target', default=None,
            help="target genome reads")
parser.add_argument('-b', '--boost', default=False,type=bool,
            help="True if new boosting scores are to be generated")
parser.add_argument('-c', '--coarse', default=11,type=int,
            help="coarse filtering kmer length")
parser.add_argument('-m', '--match', default=9,type=int,
            help="fine filtering kmer length")
parser.add_argument('-f', '--filter', default=10,type=int,
            help="Kmer mismatches allowed during fine filtering.")
parser.add_argument('-cn', '--copynumber', default=1,type=int,
            help="Copy number for training genome.")
parser.add_argument('-bo', '--boost_output', default="output/matchscores.csv",
            help="output file for matchscores")
parser.add_argument('-to', '--target_output', default="results/predictions.csv",
            help="output file for query copy number predictions")
parser.add_argument('-r', '--results', default=False,type=bool,
            help="True if target results are to be written")
opts = parser.parse_args()


# Variable declarations and imports
# # #####################################
filteringKmerLength = opts.coarse
matchingKmerLength = opts.match

f0 = opts.filter

templates = list(SeqIO.parse("templates.fasta", "fasta"))
templateNames = [t.id for t in templates]
templates = [str(t.seq) for t in templates]


#Generate k-mers
templateKmers = kmerGenerator(templates, filteringKmerLength)

if opts.boost:
    run_boosts(opts.copynumber, opts.boost_output,templates,templateNames,templateKmers,filteringKmerLength,matchingKmerLength,f0)

if opts.target != None:
    test_target([opts.target],opts.target_output,opts.boost_output, templates, templateKmers,filteringKmerLength, matchingKmerLength, f0)
