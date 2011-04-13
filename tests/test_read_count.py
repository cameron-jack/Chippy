import sys
sys.path.extend(['../src', '..'])
import warnings
warnings.filterwarnings('ignore',
        "you've sliced beyond the limits of the contig")

import numpy

from cogent import LoadTable
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files
from chippy.core.read_count import WholeChrom
from chippy.util.definition import PLUS_STRAND, MINUS_STRAND, NULL_STRAND

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
        
        for start, end, expect in coords_expect:
            got = read_counter[start:end]
            self.assertEqual(got, expect)
        
        remove_files(['sample.txt.gz'])
    
    def test_slicing_matches_numpy(self):
        """slicing a WholeChrom instance should work for +/- strands"""
        read_counter = WholeChrom('sample.txt.gz', max_read_length=30)
        got_plus = read_counter[5:20]
        expect_plus = numpy.array([5, 5, 5, 5, 5, 5, 5, 2, 0, 0, 10, 10, 10,
                                    10, 10])
        got_plus = read_counter[5:20:1]
        expect_plus = numpy.array([5, 5, 5, 5, 5, 5, 5, 2, 0, 0, 10, 10, 10,
                                    10, 10])
        self.assertEqual(got_plus, expect_plus)
        got_minus = read_counter[20:5:-1]
        expect_minus = numpy.array([10, 10, 10, 10, 10, 10, 0, 0, 2, 5, 5, 5,
                                    5, 5, 5])
        self.assertEqual(got_minus, expect_minus)
        got = read_counter[15:30]
        expect = numpy.array([10]*10)
        self.assertEqual(got, expect)
        remove_files(['sample.txt.gz'])
    
    def test_total_chrom(self):
        """whole `chromosome' should be correct"""
        read_counter = WholeChrom('sample.txt.gz', max_read_length=30)
        chrom = numpy.array([0, 0, 2, 5, 5, 5, 5, 5, 5, 5, 5, 5, 2, 0, 0,
                10, 10, 10, 10, 10, 10, 10, 10, 10, 10], numpy.int32)
        self.assertEqual(read_counter.counts, chrom)
        remove_files(['sample.txt.gz'])
    
    def test_count_reads_plus_strand(self):
        """correctly count reads from plus strand"""
        read_counter = WholeChrom('sample.txt.gz', strand=1, max_read_length=30)
        coords_expect = [(10, 15, numpy.array([3,3,0,0,0], numpy.uint32)),
                         (0, 5, numpy.array([0,0,2,3,3], numpy.uint32)),
                         (12, 16, numpy.array([0,0,0,10], numpy.uint32)),
                         (21, 25, numpy.array([10,10,10,10], numpy.uint32))]
        for start, end, expect in coords_expect:
            got = read_counter[start:end:1]
            self.assertEqual(got, expect)
        
        remove_files(['sample.txt.gz'])
    
    def test_add(self):
        """should correctly handle the += idiom"""
        read_counter = WholeChrom('sample.txt.gz', strand=-1, max_read_length=30)
        read_counter[0:10] += 4
        expect = numpy.zeros(read_counter.counts.shape[0], int)
        expect[3:13] += 2
        expect[:10] += 4
        self.assertEqual(read_counter.counts, expect)
        remove_files(['sample.txt.gz'])
    
    def test_add_two_contigs(self):
        """add two contigs together returns correct counts"""
        chrom_1_1 = WholeChrom(counts=numpy.array([0, 0, 2, 2, 2, 2, 2, 2, 2,
                        2, 2, 2]))
        
        chrom_1_2 = WholeChrom(counts=numpy.array([0, 0, 2, 2, 2, 2, 2, 2, 2,
                        2, 2, 2, 3, 3]))
        
        chrom_1_1_and_1_2 = chrom_1_1 + chrom_1_2
        self.assertEqual(chrom_1_1_and_1_2.counts,
            numpy.array([0, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3]))
        
    
    def test_adding_contigs_strand(self):
        """adding contigs preserves alike strands, otherwise should have NULL_STRAND"""
        plus1 = WholeChrom(counts=numpy.array([0, 0, 2, 2, 2, 2, 2, 2, 2,
                        2, 2, 2]), strand=PLUS_STRAND)
        plus2 = WholeChrom(counts=numpy.array([0, 0, 2, 2, 2, 2, 2, 2, 2,
                        2, 2, 2]), strand=PLUS_STRAND)
        minus1 = WholeChrom(counts=numpy.array([0, 0, 2, 2, 2, 2, 2, 2, 2,
                        2, 2, 2, 3, 3]), strand=MINUS_STRAND)
        minus2 = WholeChrom(counts=numpy.array([0, 0, 2, 2, 2, 2, 2, 2, 2,
                        2, 2, 2]), strand=MINUS_STRAND)
        
        plusplus = plus1 + plus2
        self.assertEqual(plusplus.strand, PLUS_STRAND)
        minusminus = minus1 + minus2
        self.assertEqual(minusminus.strand, MINUS_STRAND)
        plusminus = plus1 + minus1
        self.assertEqual(plusminus.strand, NULL_STRAND)
        minusplus = minus1 + plus1
        self.assertEqual(minusplus.strand, NULL_STRAND)
    
    def test_count_reads_minus_strand(self):
        """correctly count reads from minus strand"""
        read_counter = WholeChrom(mapped_read_path='sample.txt.gz', strand=-1,
                                    max_read_length=30)
        coords_expect = [(10, 15, numpy.array([2,2,2,0,0], numpy.uint32)),
                         (0, 5, numpy.array([0,0,0,2,2], numpy.uint32)),
                         (12, 16, numpy.array([2,0,0,0], numpy.uint32)),
                         (21, 25, numpy.array([0,0,0,0], numpy.uint32))]
        for start, end, expect in coords_expect:
            got = read_counter[start:end]
            self.assertEqual(got, expect)
        
        remove_files(['sample.txt.gz'])
    

if __name__ == "__main__":
    main()
