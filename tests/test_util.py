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

if __name__ == "__main__":
    main()