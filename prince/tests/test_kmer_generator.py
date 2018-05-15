import unittest

from prince.kmer_generator import kmer_generator

class KmerGeneratorTest(unittest.TestCase):
    
    def test_kmer_generator(self):
        self.assertEqual(kmer_generator("", 0), {})
        self.assertEqual(kmer_generator("A",1), {0: ['A']})
