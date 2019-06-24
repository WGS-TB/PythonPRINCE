import unittest

import json
from Bio import SeqIO
from prince.boost import run_boosts
from pkg_resources import resource_filename
from prince.kmer_generator import kmer_generator

PATH_PREFIX = "prince/tests/"

class opts_object:
            def __init__(self, boosting_file, num_procs=2, k=25, q=20):
                self.num_procs=num_procs
                self.boosting_file = boosting_file
                self.k = k
                self.q = q

class BoostingTest(unittest.TestCase):
    
    def setUp(self):
        self.PATH_PREFIX = PATH_PREFIX
        
        self.opts1 = opts_object(PATH_PREFIX + "training_test.txt", 2)
        
        templates = list(SeqIO.parse(resource_filename('prince.resources', "templates.fasta"), "fasta"))
        templateNames = [t.id for t in templates]
        templates = [str(t.seq) for t in templates]

        #Generate k-mers
        templateKmers = kmer_generator(templates, 25)
        self.template_obj = {"Names":templateNames, "Sequences":templates, "Kmers":templateKmers}
        
        with open(resource_filename('prince.resources', "TB_primers_extended.json")) as primers:
            self.primers=json.load(primers)
        
        
    def test_run_boosts(self):
        boost_output1 = run_boosts(self.opts1, self.template_obj, self.primers)
        output_valid = True if "prince/tests/data/small_test,154,0,1,0.0" in boost_output1 else False
        assert(output_valid)
        output_valid = True if "prince/tests/data/zip_test,2163b,15,8,0.444444444444" in boost_output1 else False
        assert(output_valid)
        output_invalid = True if "random cat" in boost_output1 else False
        assert(not output_invalid)
        
        
if __name__ == '__main__':
    unittest.main()