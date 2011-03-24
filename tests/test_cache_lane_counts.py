import sys
sys.path.append('../src')

import numpy

from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files
from chippy.core.read_count import CacheLaneCounts, RegularRegionCounts

def generate_test_files():
    for index in range(10):
        counter = RegularRegionCounts(10, one_based=False)
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
        """correctly return cached counts"""
        generate_test_files()
        count_cache = CacheLaneCounts(7, '.')
        for index in range(10):
            counter = RegularRegionCounts(10, one_based=False)
            counter.addRead(index, index+1)
            self.assertEqual(count_cache[index+1],
                             counter._get_subregion_counts([(5,1)], 5))
        remove_files(count_cache.count_filenames)
    
    def test_return_centred_counts(self):
        """correctly return cached counts centred"""
        length = 10
        counter = RegularRegionCounts(length, one_based=False)
        for i in range(length):
            counter.addRead(i, length)
        
        counter.save('s_7_chr7', [(5, 1)], 5)
        
        window_size = 3
        count_cache = CacheLaneCounts(7, '.', window_size=window_size)
        self.assertEqual(count_cache[7].shape[1], 2*window_size)
        self.assertEqual(count_cache[7][0], range(window_size, 3*window_size))
        remove_files(count_cache.count_filenames)
    

if __name__ == "__main__":
    main()