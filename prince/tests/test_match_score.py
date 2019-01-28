import unittest

from prince.match_score import get_reads_records

class OpenFiles(unittest.TestCase):
    
    def setUp(self):
        self.PATH_PREFIX = "prince/tests/data/"
        
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

if __name__ == '__main__':
    unittest.main()