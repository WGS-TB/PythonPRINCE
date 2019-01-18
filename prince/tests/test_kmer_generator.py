import unittest

from prince.kmer_generator import kmer_generator

class KmerGeneratorTest(unittest.TestCase):
    
    def test_kmer_generator(self):
        self.assertEqual(kmer_generator("", 0), {})
        self.assertEqual(kmer_generator("A",1), {0: ['A', 'T']})
        self.assertEqual(kmer_generator(["ATTCG","GGGG"],3), {0: ['ATT', 'TTC', 'TCG', 'CGA', 'GAT', 'ATC', 'TCG', 'CGA', 'GAA', 'AAT'], 1: ['GGG', 'GGG', 'GGG', 'GGG', 'CCC', 'CCC', 'CCC', 'CCC']})
        self.assertEqual(kmer_generator(["ATTCG","GGGG"],3, extension=False), {0: ['ATT', 'TTC', 'TCG', 'CGA', 'GAA', 'AAT'], 1: ['GGG', 'GGG', 'CCC', 'CCC']})
