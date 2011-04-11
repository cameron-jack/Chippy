import sys
sys.path.extend(['../src', '..'])

import numpy

from cogent import LoadTable
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files
from chippy.core.read_count import WholeChrom

# sample counts for below
# Read num / pos    0   1   2   3   4   5   6   7   8   9   10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25
#                 0 .   .   2   2   2   2   2   2   2   2   2   2   .   .   .   .   .   .   .   .   .   .   .   .   .   .
#                 1 .   .       1   1   1   1   1   1   1   1   1   .   .   .   .   .   .   .   .   .   .   .   .   .   .
#     (- strand)  2 .   .   .   2   2   2   2   2   2   2   2   2   2   .   .   .   .   .   .   .   .   .   .   .   .   .
#                 3 .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   10  10  10  10  10  10  10  10  10  10  .
#         counts    0   0   2   5   5   5   5   5   5   5   5   5   2   0   0   10  10  10  10  10  10  10  10  10  10  0


class MinimalRegionCountTests(TestCase):
    def setUp(self):
        """write some read data"""
        header = ['start', 'length', 'strand', 'freq']
        mapped = [[2, 10, 1, 2],
                  [3, 9, 1, 1],
                  [3, 10, -1, 2],
                  [15, 10, 1, 10]]
        read_table = LoadTable(header=header, rows=mapped)
        read_table.writeToFile('sample.txt.gz', sep='\t')
    
    def test_count_reads(self):
        """correctly count reads from both strands"""
        read_counter = WholeChrom('sample.txt.gz', max_read_length=30)
        coords_expect = [(10, 15, numpy.array([5,5,2,0,0], numpy.uint32)),
                         (0, 5, numpy.array([0,0,2,5,5], numpy.uint32)),
                         (12, 16, numpy.array([2,0,0,10], numpy.uint32)),
                 (21, 25, numpy.array([10,10,10,10], numpy.uint32))]
        
        read_counter.update()
        for start, end, expect in coords_expect:
            got = read_counter[start:end]
            self.assertEqual(got, expect)
        
        remove_files(['sample.txt.gz'])
    
    def test_total_chrom(self):
        """whole `chromosome' should be correct"""
        read_counter = WholeChrom('sample.txt.gz', max_read_length=30)
        read_counter.update()
        chrom = numpy.array([0, 0, 2, 5, 5, 5, 5, 5, 5, 5, 5, 5, 2, 0, 0,
                10, 10, 10, 10, 10, 10, 10, 10, 10, 10], numpy.int32)
        self.assertEqual(read_counter.total_count, chrom)
    
    def test_count_reads_plus_strand(self):
        """correctly count reads from plus strand"""
        read_counter = WholeChrom('sample.txt.gz', strand=1, max_read_length=30)
        coords_expect = [(10, 15, numpy.array([3,3,0,0,0], numpy.uint32)),
                         (0, 5, numpy.array([0,0,2,3,3], numpy.uint32)),
                         (12, 16, numpy.array([0,0,0,10], numpy.uint32)),
                         (21, 25, numpy.array([10,10,10,10], numpy.uint32))]
        read_counter.update()
        for start, end, expect in coords_expect:
            got = read_counter[start:end:1]
            self.assertEqual(got, expect)
        
        remove_files(['sample.txt.gz'])
    
    def test_add(self):
        """should correctly handle the += idiom"""
        read_counter = WholeChrom('sample.txt.gz', strand=-1, max_read_length=30)
        read_counter[0:10] += 4
        expect = numpy.zeros(read_counter.total_count.shape[0], int)
        expect[:10] += 4
        self.assertEqual(read_counter.total_count, expect)
        remove_files(['sample.txt.gz'])
    
    def test_count_reads_minus_strand(self):
        """correctly count reads from minus strand"""
        read_counter = WholeChrom('sample.txt.gz', strand=-1, max_read_length=30)
        coords_expect = [(10, 15, numpy.array([2,2,2,0,0], numpy.uint32)),
                         (0, 5, numpy.array([0,0,0,2,2], numpy.uint32)),
                         (12, 16, numpy.array([2,0,0,0], numpy.uint32)),
                         (21, 25, numpy.array([0,0,0,0], numpy.uint32))]
        
        read_counter.update()
        for start, end, expect in coords_expect:
            got = read_counter[start:end]
            self.assertEqual(got, expect)
        
        remove_files(['sample.txt.gz'])
    
    def test_compare_old_new(self):
        """RegularRegionCounts and WholeChrom should give same answer"""
        read_counter = WholeChrom('sample.txt.gz', max_read_length=30)
        read_counter.update()
        counter = RegularRegionCounts(25, one_based=False)
        
        for i in range(2):
            counter.addRead(2, 12)
            counter.addRead(3, 13)
        counter.addRead(3, 12)
        
        for i in range(10):
            counter.addRead(15, 25)
        
        self.assertEqual(counter.counts, read_counter.total_count)
    

if __name__ == "__main__":
    main()
