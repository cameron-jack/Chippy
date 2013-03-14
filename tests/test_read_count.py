import sys
sys.path.extend(['../src', '..'])
import warnings
warnings.filterwarnings('ignore',
        "Not using MPI as mpi4py not found")

from chippy.core.region_of_interest import ROI
from cogent.util.unit_test import TestCase
from chippy.core.count_tags import ROI
from chippy.core.read_count import read_BAM, read_BED
from chippy.express.db_schema import Gene

from chippy.util.definition import PLUS_STRAND, MINUS_STRAND, NULL_STRAND

class MinimalRegionCountTests(TestCase):
    """
    tests for: add_counts_to_ROI, read_BAM, read_BED,ROI.
    """

    def setUpROIsForSynthTests(self):
        """
        Minimal regions that should be easy to test
        EVERYTHING is in 0-based space.
        """
        # Actual gene end=10 but in Python 0-1 [,) space we add one so we can slice it
        g1 = Gene('g1', 'g1', 'protein_coding', 'fake_g1', 'ok', '1', 4, 11, PLUS_STRAND)
        g2 = Gene('g2', 'g2', 'protein_coding', 'fake_g2', 'ok', '1', 4, 11, MINUS_STRAND)

        window_radius = 5
        ROIs = []

        # window should be indices [-1,0,1,2,3, 4,5,6,7,8]
        win_start, win_end = g1.getTssCentredCoords(window_radius)
        roi = ROI(g1, win_start, win_end)
        ROIs.append(roi)

        # window should be indices [15,14,13,12,11, 10,9,8,7,6]
        win_start, win_end = g2.getTssCentredCoords(window_radius)
        roi = ROI(g2, win_start, win_end)
        ROIs.append(roi)

        window_radius=10
        # window should be indices [20,19,18,17,16,15,14,13,12,11, 10,9,8,7,6,5,4,3,2,1]
        win_start, win_end = g2.getTssCentredCoords(window_radius)
        roi = ROI(g2, win_start, win_end)
        ROIs.append(roi)

        return ROIs

    def test_add_counts_to_ROI(self):
        """ Source irrespective, reads should map onto ROI correcly """
        ROIs = self.setUpROIsForSynthTests()

        # length 1 reads with positive strand
        entry_start=4 # At TSS
        entry_end=5
        ROIs[0].add_counts_to_ROI(entry_start, entry_end)
        self.assertEqual(ROIs[0].counts, [0, 0, 0, 0, 0, 1, 0, 0, 0, 0])
        entry_start=8
        entry_end=9
        ROIs[0].add_counts_to_ROI(entry_start, entry_end)
        self.assertEqual(ROIs[0].counts, [0, 0, 0, 0, 0, 1, 0, 0, 0, 1])
        entry_start=9 # this is to the right of the ROI
        entry_end=10
        ROIs[0].add_counts_to_ROI(entry_start, entry_end)
        self.assertEqual(ROIs[0].counts, [0, 0, 0, 0, 0, 1, 0, 0, 0, 1])
        entry_start=-1 # should work okay with negative indices
        entry_end=0
        ROIs[0].add_counts_to_ROI(entry_start, entry_end)
        self.assertEqual(ROIs[0].counts, [1, 0, 0, 0, 0, 1, 0, 0, 0, 1])
        entry_start=-2 # this is to the left of the ROI
        entry_end=-1
        ROIs[0].add_counts_to_ROI(entry_start, entry_end)
        self.assertEqual(ROIs[0].counts, [1, 0, 0, 0, 0, 1, 0, 0, 0, 1])

        # length 1 reads with negative strand
        entry_start=10 # At TSS
        entry_end=11
        ROIs[1].add_counts_to_ROI(entry_start, entry_end)
        self.assertEqual(ROIs[1].counts, [0, 0, 0, 0, 0, 1, 0, 0, 0, 0])
        entry_start=6 # RHS of window
        entry_end=7
        ROIs[1].add_counts_to_ROI(entry_start, entry_end)
        self.assertEqual(ROIs[1].counts, [0, 0, 0, 0, 0, 1, 0, 0, 0, 1])
        entry_start=5 # 1 beyond RHS of window
        entry_end=6
        ROIs[1].add_counts_to_ROI(entry_start, entry_end)
        self.assertEqual(ROIs[1].counts, [0, 0, 0, 0, 0, 1, 0, 0, 0, 1])
        entry_start=15 # LHS of window
        entry_end=16
        ROIs[1].add_counts_to_ROI(entry_start, entry_end)
        self.assertEqual(ROIs[1].counts, [1, 0, 0, 0, 0, 1, 0, 0, 0, 1])
        entry_start=16 # 1 beyond LHS of window
        entry_end=17
        ROIs[1].add_counts_to_ROI(entry_start, entry_end)
        self.assertEqual(ROIs[1].counts, [1, 0, 0, 0, 0, 1, 0, 0, 0, 1])


        # First case, perfect overlap, slightly inside the ROI
        entry_start = 1
        entry_end = 8 # end and include genome position 7
        ROIs[0].add_counts_to_ROI(entry_start, entry_end)
        self.assertEqual(ROIs[0].counts, [1, 0, 1, 1, 1, 2, 1, 1, 1, 1])
        ROIs[1].add_counts_to_ROI(entry_start, entry_end)
        self.assertEqual(ROIs[1].counts, [1, 0, 0, 0, 0, 1, 0, 0, 1, 2])
        ROIs[2].add_counts_to_ROI(entry_start, entry_end)
        self.assertEqual(ROIs[2].counts, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 1, 1, 1, 1, 1, 1, 1])

        # test both addition to previous ROI and positive over-the-edge counts
        entry_start = 7
        entry_end = 13
        ROIs[0].add_counts_to_ROI(entry_start, entry_end)
        self.assertEqual(ROIs[0].counts, [1, 0, 1, 1, 1, 2, 1, 1, 2, 2])
        ROIs[1].add_counts_to_ROI(entry_start, entry_end)
        self.assertEqual(ROIs[1].counts, [1, 0, 0, 1, 1, 2, 1, 1, 2, 2])
        ROIs[2].add_counts_to_ROI(entry_start, entry_end)
        self.assertEqual(ROIs[2].counts, [0, 0, 0, 0, 0, 0, 0, 0, 1, 1,  1, 1, 1, 2, 1, 1, 1, 1, 1, 1])

        # test addition and negative over-the-edge counts
        entry_start = -2
        entry_end = 5
        ROIs[0].add_counts_to_ROI(entry_start, entry_end)
        self.assertEqual(ROIs[0].counts, [2, 1, 2, 2, 2, 3, 1, 1, 2, 2])
        ROIs[1].add_counts_to_ROI(entry_start, entry_end)
        self.assertEqual(ROIs[1].counts, [1, 0, 0, 1, 1, 2, 1, 1, 2, 2])
        ROIs[2].add_counts_to_ROI(entry_start, entry_end)
        self.assertEqual(ROIs[2].counts, [0, 0, 0, 0, 0, 0, 0, 0, 1, 1,  1, 1, 1, 2, 1, 1, 2, 2, 2, 2])

    def setUpROIsForFiles(self):
        """
        In brca2-11.sam if we make TSS= 151341224 and window radius = 5.
        Then we should get 119->223 at 0 and 224->229 at +2.
        If we take TSS= 151345949 and window radius = 5, then we should get
        945->949 at +2 and 950->954 at 0, in 1-offset space.
        """
        g1 = Gene('BRCA2a', 'brca2a', 'protein_coding', 'fake_brca2',
                'ok', 'chr_5', 151341223, 151345949, PLUS_STRAND)
        g2 = Gene('BRCA2b', 'brca2b', 'protein_coding', 'fake_brca2',
                'ok', 'chr_5', 151341223, 151345949, MINUS_STRAND)
        g3 = Gene('BRCA2c', 'brca2c', 'protein_coding', 'fake_brca2',
                'ok', 'chr_11', 151341223, 151345949, PLUS_STRAND)

        ROIs = []

        # [0, 0, 0, 0, 0, 2, 2, 2, 2, 2]
        window_radius = 5
        win_start, win_end = g1.getTssCentredCoords(window_radius)
        roi1 = ROI(g1, win_start, win_end)
        ROIs.append(roi1)

        window_radius = 5
        # [0, 0, 0, 0, 0, 0, 2, 2, 2, 2] due to reads finishing 1 early
        win_start, win_end = g2.getTssCentredCoords(window_radius)
        roi2 = ROI(g2, win_start, win_end)
        ROIs.append(roi2)

        # perfectly overlapping, and larger than roi1
        # This tests that the readers can handle overlapping ROIs
        # It also shows that the window shifts properly
        # [0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2]
        window_radius = 6
        win_start, win_end = g1.getTssCentredCoords(window_radius)
        roi3 = ROI(g1, win_start, win_end)
        ROIs.append(roi3)

        # Minus strand version of above
        # [0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2]
        window_radius = 6
        win_start, win_end = g2.getTssCentredCoords(window_radius)
        roi4 = ROI(g2, win_start, win_end)
        ROIs.append(roi4)

        # Test alternative chromosomes of interest
        # [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        window_radius = 6
        win_start, win_end = g3.getTssCentredCoords(window_radius)
        roi5 = ROI(g3, win_start, win_end)
        ROIs.append(roi5)

        # NOTE: No alternative chroms in test BED file to test for...

        return ROIs

    def test_count_BAM(self):
        """ make sure we get the expected counts from BAM reading.
        Negative stand expectations are offset here because that's what's
        in the data """
        ROIs = self.setUpROIsForFiles()
        ROIs, num_tags, num_bases, rr = read_BAM('data/brca2-11_sorted.bam',
                ROIs)

        self.assertEqual(num_tags, 8) # 2 reads for each of the first 4 cases
        self.assertEqual(num_bases, 40) # sum of all counts
        self.assertEqual(ROIs[0].counts, [0, 0, 0, 0, 0, 2, 2, 2, 2, 2])
        self.assertEqual(ROIs[1].counts, [0, 0, 0, 0, 0, 0, 2, 2, 2, 2])
        self.assertEqual(ROIs[2].counts, [0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2])
        self.assertEqual(ROIs[3].counts, [0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2])
        self.assertEqual(ROIs[4].counts, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

    def test_count_BED(self):
        """ make sure we get the expected counts from BED reading
        read_BED will often mess up the order of the ROIs so we'll need
        a clever method to index them to check them correctly.
        Negative stand expectations are offset here because that's what's
        in the data. """
        ROIs = self.setUpROIsForFiles()
        ROIs, num_tags, num_bases, rr = read_BED('data/brca2-11.bed', ROIs)
        self.assertEqual(num_tags, 8) # 2 reads for each of the first 4 cases
        self.assertEqual(num_bases, 40) # sum of all counts

        # read_BED will often mess up the order of the ROIs so we'll need
        # a clever method to index them to check them correctly

        for roi in ROIs:
            if len(roi.counts) == 10 and sum(roi.counts) == 10:
                self.assertEqual(roi.counts, [0, 0, 0, 0, 0, 2, 2, 2, 2, 2])
            elif len(roi.counts) == 10 and sum(roi.counts) == 8:
                self.assertEqual(roi.counts, [0, 0, 0, 0, 0, 0, 2, 2, 2, 2])
            elif len(roi.counts) == 12 and sum(roi.counts) == 12:
                self.assertEqual(roi.counts, [0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2])
            elif len(roi.counts) == 12 and sum(roi.counts) == 10:
                self.assertEqual(roi.counts, [0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2])
            elif len(roi.counts) == 12 and sum(roi.counts) == 0:
                self.assertEqual(roi.counts, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
            else:
                self.assertEqual(roi.counts, roi.gene_id)
