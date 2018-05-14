import unittest

from prince.kmer_generator import kmerGenerator

class KmerGeneratorTest(unittest.TestCase):
    
    def testKmerGenerator(self):
        self.assertEqual(kmerGenerator("", 0), {})
        self.assertEqual(kmerGenerator("A",1), {0: ['A']})
