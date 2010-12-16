import sys
sys.path.append('../src')

import util

import numpy

from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files, flatten


class GroupingTests(TestCase):
    def test_no_remainder(self):
        """return the correct groups when there's no remainder"""
        data = range(20)
        grouped = util.make_even_groups(data, 5)
        self.assertEqual(len(grouped), 4)
        for group in grouped:
            self.assertEqual(len(group), 5)
        full = sorted(flatten(grouped))
        self.assertEqual(full, data)
    
    def test_with_remainder(self):
        """return the correct groups when there's a remainder"""
        data = range(21)
        grouped = util.make_even_groups(data, 5)
        self.assertEqual(len(grouped), 4)
        for group in grouped:
            self.assertEqual(len(group), 5)
        full = sorted(flatten(grouped))
        self.assertEqual(full, data[:-1])
    
    def test_one_group(self):
        """only one group"""
        data = range(20)
        grouped = util.make_even_groups(data, 20)
        self.assertEqual(len(grouped), 1)
        for group in grouped:
            self.assertEqual(len(group), 20)
        full = sorted(flatten(grouped))
        self.assertEqual(full, data)
    
    def test_centred_coords(self):
        """should correctly return the central fragment"""
        start, end = util.get_centred_coords(100, 10)
        self.assertEqual(start, 40)
        self.assertEqual(end, 60)
        self.assertEqual(util.get_centred_coords(100, 50), (0, 100))
        self.assertRaises(AssertionError, util.get_centred_coords, 100, 60)
    
    def test_unique_records(self):
        """unique_records returns correctly ordered unique data """
        data = range(20)
        for i in (11, 15, 18):
            data.insert(i, 10) # inserts redundant elements
        got = util.unique_records(data)
        self.assertEqual(got, range(20))
        

if __name__ == "__main__":
    main()