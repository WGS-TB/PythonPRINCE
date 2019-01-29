import unittest

from prince.predict import get_copy_number

class GetCopyNumberTest(unittest.TestCase):
    
    def test_get_copy_number(self):
        self.assertEqual(get_copy_number(0, 0, 0), 0)
        self.assertEqual(get_copy_number(1, 0, 0), 0)
        self.assertEqual(get_copy_number(0, 1, 0), 0)
        self.assertEqual(get_copy_number(0, 0, 1), 1)
        self.assertEqual(get_copy_number(0, 1, 1), 1)
        self.assertEqual(get_copy_number(1, 0, 1), 1)
        self.assertEqual(get_copy_number(1, 1, 0), 1)
        self.assertEqual(get_copy_number(1, 1, 1), 2)
        self.assertEqual(get_copy_number(-1, 0, 0), 0)
        self.assertEqual(get_copy_number(0, -1, 0), 0)
        self.assertEqual(get_copy_number(0, 0, -1), -1)
        self.assertEqual(get_copy_number(0, -1, -1), -1)
        self.assertEqual(get_copy_number(-1, 0, -1), -1)
        self.assertEqual(get_copy_number(-1, -1, 0), 1)
        self.assertEqual(get_copy_number(-1, -1, -1), 0)
        self.assertEqual(get_copy_number(0, 1, -1), -1)
        self.assertEqual(get_copy_number(0, -1, 1), 1)
        self.assertEqual(get_copy_number(-1, 0, 1), 1)
        self.assertEqual(get_copy_number(100, 100, 100), 10100)
        self.assertEqual(get_copy_number(987, 896, 876), 885228)
        
        
if __name__ == '__main__':
    unittest.main()