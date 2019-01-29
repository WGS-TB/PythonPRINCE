import unittest

import json
from Bio import SeqIO
from prince.match_score import get_reads_records, compute_match_score
from pkg_resources import resource_filename
from prince.kmer_generator import kmer_generator

PATH_PREFIX = "prince/tests/data/"


class OpenFiles(unittest.TestCase):

    def setUp(self):
        self.PATH_PREFIX = PATH_PREFIX 

    def test_nonsense_record(self):
        with self.assertRaises(IOError):
            get_reads_records("dasd")
    
    def test_get_small_reads_records(self):
        record1,record2,gzip_handle1,gzip_handle2 = get_reads_records(self.PATH_PREFIX + "small_test")
        self.assertIsNotNone(record1)
        self.assertIsNotNone(record2)
        self.assertNotEqual(record1, record2)
        self.assertIsNone(gzip_handle1)
        self.assertIsNone(gzip_handle2)
        
    def test_get_medium_reads_records(self):
        record1,record2,gzip_handle1,gzip_handle2 = get_reads_records(self.PATH_PREFIX + "medium_test")
        self.assertIsNotNone(record1)
        self.assertIsNotNone(record2)
        self.assertNotEqual(record1, record2)
        self.assertIsNone(gzip_handle1)
        self.assertIsNone(gzip_handle2)
        
    def test_get_zip_reads_records(self):
        record1,record2,gzip_handle1,gzip_handle2 = get_reads_records(self.PATH_PREFIX + "zip_test")
        self.assertIsNotNone(record1)
        self.assertIsNotNone(record2)
        self.assertNotEqual(record1, record2)
        self.assertNotEqual(gzip_handle1, gzip_handle2)
        self.assertIsNotNone(gzip_handle1)
        self.assertIsNotNone(gzip_handle2)
        
    def test_get_zip_two_lines_reads_records(self):
        record1,record2,gzip_handle1,gzip_handle2 = get_reads_records(self.PATH_PREFIX + "zip_test2_1.fastq.gz", self.PATH_PREFIX + "zip_test2_2.fastq.gz")
        self.assertIsNotNone(record1)
        self.assertIsNotNone(record2)
        self.assertNotEqual(record1, record2)
        self.assertNotEqual(gzip_handle1, gzip_handle2)
        self.assertIsNotNone(gzip_handle1)
        self.assertIsNotNone(gzip_handle2)
        
        
class MatchScoreTest(unittest.TestCase):
    
    def setUp(self):
        self.PATH_PREFIX = PATH_PREFIX
        
        templates = list(SeqIO.parse(resource_filename('prince.resources', "templates.fasta"), "fasta"))
        templateNames = [t.id for t in templates]
        templates = [str(t.seq) for t in templates]

        #Generate k-mers
        templateKmers = kmer_generator(templates, 25)
        self.template_obj = {"Names":templateNames, "Sequences":templates, "Kmers":templateKmers}
        
        with open(resource_filename('prince.resources', "TB_primers_extended.json")) as primers:
            self.primers=json.load(primers)
        
        
    def test_compute_match_score(self):
        match_score_small = compute_match_score(self.PATH_PREFIX + "small_test", "", self.template_obj, 25, self.primers)
        match_score_medium_1 = compute_match_score(self.PATH_PREFIX + "medium_test", "", self.template_obj, 25, self.primers)
        
        match_score_medium_2 = compute_match_score(self.PATH_PREFIX + "medium_test1.fq", self.PATH_PREFIX + "medium_test2.fq", self.template_obj, 25, self.primers)
        
        self.assertEqual(match_score_small, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        self.assertEqual(match_score_medium_1, match_score_medium_2)
        self.assertEqual(match_score_medium_1, [1.0, 3.0, 2.0, 1.0, 0.0, 1.0, 2.0, 4.0, 1.0, 3.0, 0.5, 0.0, 0.0, 7.0, 0.25, 0.4444444444444444, 0.0, 4.0, 0.0, 0.0, 0.0, 8.0, 6.0, 0.0])

if __name__ == '__main__':
    unittest.main()