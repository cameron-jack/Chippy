import sys
sys.path.append('../src')

import numpy

from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files
from region_count import CacheLaneCounts,RegionCounts

def generate_test_files():
    for index in range(10):
        counter = RegionCounts(10, one_based=False)
        counter.addRead(index, index+1)
        counter.save('s_7_chr%d'%(index+1), [(5,1)], 5)


class CacheLaneCountTests(TestCase):
    def test_proper_initialise(self):
        generate_test_files()
        try:
            count_cache = CacheLaneCounts(7, '../src')
        except IOError:
            pass
            #print 'No Files There!'

        try:
            count_cache = CacheLaneCounts(8, '.')
        except IOError:
            pass
            #print 'No Lane 8 Files There!'

        try:
            count_cache = CacheLaneCounts(8, '../dump')
        except IOError:
            pass
            #print 'Not a Valid Path!'

        count_cache = CacheLaneCounts(7, '.')
        self.assertEqual(count_cache.num_files, 10)
        self.assertEqual(count_cache.lane, 7)
        self.assertEqual(count_cache.path, '.')
        remove_files(count_cache.count_filenames)

    def test_return_counts(self):
        generate_test_files()
        count_cache = CacheLaneCounts(7, '.')
        for index in range(10):
            counter = RegionCounts(10, one_based=False)
            counter.addRead(index, index+1)
            self.assertEqual(count_cache.getCountsForChrom(index+1),
                             counter._get_subregion_counts([(5,1)], 5))
        remove_files(count_cache.count_filenames)

if __name__ == "__main__":
    main()