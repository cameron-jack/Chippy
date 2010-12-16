import sys
sys.path.append('../src')

import numpy

from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files
from region_count import RegionCounts

class RegionCountTests(TestCase):
    def test_handle_one_zero_based_systems(self):
        """should correctly handle zero or one based counting system"""
        count_start = 1
        count_end = 10
        counter = RegionCounts(100, one_based=True)
        counter.addRead(count_start, count_end)
        self.assertEqual(counter.counts[0:10], numpy.ones(10, numpy.uint16))
        
        counter = RegionCounts(100, one_based=False)
        counter.addRead(count_start, count_end)
        expect = numpy.ones(10, numpy.uint16)
        expect[0] = 0
        self.assertEqual(counter.counts[0:10], expect)
    
    def test_subregion_counts(self):
        """should return correct counts fo subregions"""
        counter = RegionCounts(100, one_based=False)
        
        counter.addRead(20, 30)
        counter.addRead(75, 85)
        
        tss_strand = [(20, 1), (84, -1)]
        
        subregions = counter._get_subregion_counts(tss_strand, 10)
        self.assertEqual(subregions.shape, (2, 20))
        col_sums = subregions.sum(axis=0)
        expect = numpy.zeros(20, numpy.uint16)
        expect[10:] = 2
        self.assertEqual(col_sums, expect)
    
    def test_save_counts(self):
        """correctly saves to a file"""
        file_name = 'saved_counts'
        counter = RegionCounts(100, one_based=False)
        counter.addRead(20, 30)
        counter.addRead(75, 85)
        tss_strand = [(20, 1), (84, -1)]
        expect = counter._get_subregion_counts(tss_strand, 10)
        counter.save(file_name, tss_strand, 10)
        got = numpy.load(file_name+'.npy')
        self.assertEqual(got, expect)
        remove_files([file_name+'.npy'])
    

if __name__ == "__main__":
    main()