import sys
sys.path.extend(['../src', '..'])
import warnings
warnings.filterwarnings('ignore',
        "Not using MPI as mpi4py not found")

from cogent.util.unit_test import TestCase, main
from chippy.core.count_tags import ROI
from chippy.core.read_count import read_BAM, read_BED, add_counts_to_ROI
from chippy.express.db_schema import Gene

from chippy.util.definition import PLUS_STRAND, MINUS_STRAND, NULL_STRAND

class MinimalRegionCountTests(TestCase):
    """
    tests for: add_counts_to_ROI, read_BAM, read_BED,ROI.
    """

    def setUpROIsForSynthTests(self):
        """
        Minimal regions that should be easy to test
        TSS is in 1-based space
        """
        g1 = Gene('g1', 'g1', 'protein_coding', 'fake_g1', 'ok', '1', 4, 10, PLUS_STRAND)
        g2 = Gene('g2', 'g2', 'protein_coding', 'fake_g2', 'ok', '1', 4, 10, MINUS_STRAND)

        window_size = 5
        ROIs = []

        # window should be indices [_,_,1,2,3, 4,5,6,7,8]
        win_start, win_end = g1.getTssCentredCoords(window_size)
        roi = ROI(g1, win_start, win_end)
        ROIs.append(roi)

        # window should be indices [15,14,13,12,11, 10,9,8,7,6]
        win_start, win_end = g2.getTssCentredCoords(window_size)
        roi = ROI(g2, win_start, win_end)
        ROIs.append(roi)

        window_size=10
        # window should be indices [20,19,18,17,16,15,14,13,12,11, 10,9,8,7,6,5,4,3,2,1]
        win_start, win_end = g2.getTssCentredCoords(window_size)
        roi = ROI(g2, win_start, win_end)
        ROIs.append(roi)

        return ROIs

    def test_add_counts_to_ROI(self):
        """ If this function works, the rest should easy """
        ROIs = self.setUpROIsForSynthTests()

        # First case, perfect overlap, slightly inside the ROI
        entry_start = 1 # start at genome position 1
        entry_end = 7 # end and include genome position 7
        add_counts_to_ROI(ROIs[0], entry_start, entry_end)
        self.assertEqual(ROIs[0].counts, [0, 0, 1, 1, 1, 1, 1, 1, 1, 0])
        add_counts_to_ROI(ROIs[1], entry_start, entry_end)
        self.assertEqual(ROIs[1].counts, [0, 0, 0, 0, 0, 0, 0, 0, 1, 1])
        add_counts_to_ROI(ROIs[2], entry_start, entry_end)
        self.assertEqual(ROIs[2].counts, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 1, 1, 1, 1, 1, 1, 1])

        # test both addition to previous ROI and positive over-the-edge counts
        entry_start = 7
        entry_end = 12
        add_counts_to_ROI(ROIs[0], entry_start, entry_end)
        self.assertEqual(ROIs[0].counts, [0, 0, 1, 1, 1, 1, 1, 1, 2, 1])
        add_counts_to_ROI(ROIs[1], entry_start, entry_end)
        self.assertEqual(ROIs[1].counts, [0, 0, 0, 1, 1, 1, 1, 1, 2, 1])
        add_counts_to_ROI(ROIs[2], entry_start, entry_end)
        self.assertEqual(ROIs[2].counts, [0, 0, 0, 0, 0, 0, 0, 0, 1, 1,  1, 1, 1, 2, 1, 1, 1, 1, 1, 1])

        # test addition and negative over-the-edge counts
        entry_start = -2
        entry_end = 2
        add_counts_to_ROI(ROIs[0], entry_start, entry_end)
        self.assertEqual(ROIs[0].counts, [1, 1, 2, 2, 1, 1, 1, 1, 2, 1])
        add_counts_to_ROI(ROIs[1], entry_start, entry_end)
        self.assertEqual(ROIs[1].counts, [0, 0, 0, 1, 1, 1, 1, 1, 2, 1])
        add_counts_to_ROI(ROIs[2], entry_start, entry_end)
        self.assertEqual(ROIs[2].counts, [0, 0, 0, 0, 0, 0, 0, 0, 1, 1,  1, 1, 1, 2, 1, 1, 1, 1, 2, 2])

    def setUpROIsForFiles(self):
        """
        In brca2-11.sam if we make TSS= 151341224 and window size = 5.
        Then we should get 119->223 at 0 and 224->229 at +2.
        If we take TSS= 151345949 and window size = 5, then we should get
        945->949 at +2 and 950->954 at 0, in 1-offset space.
        """
        g1 = Gene('BRCA2a', 'brca2a', 'protein_coding', 'fake_brca2', 'ok', 'chr_5', 151341224, 151345949, PLUS_STRAND)
        g2 = Gene('BRCA2b', 'brca2b', 'protein_coding', 'fake_brca2', 'ok', 'chr_5', 151341224, 151345949, MINUS_STRAND)
        g3 = Gene('BRCA2c', 'brca2c', 'protein_coding', 'fake_brca2', 'ok', 'chr_11', 151341224, 151345949, PLUS_STRAND)

        ROIs = []

        # [0, 0, 0, 0, 0, 2, 2, 2, 2, 2]
        window_size = 5
        win_start, win_end = g1.getTssCentredCoords(window_size)
        roi1 = ROI(g1, win_start, win_end)
        ROIs.append(roi1)

        window_size = 5
        # [0, 0, 0, 0, 0, 0, 2, 2, 2, 2] due to reads finishing 1 early
        win_start, win_end = g2.getTssCentredCoords(window_size)
        roi2 = ROI(g2, win_start, win_end)
        ROIs.append(roi2)

        # perfectly overlapping, and larger than roi1
        # This tests that the readers can handle overlapping ROIs
        # It also shows that the window shifts properly
        # [0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2]
        window_size = 6
        win_start, win_end = g1.getTssCentredCoords(window_size)
        roi3 = ROI(g1, win_start, win_end)
        ROIs.append(roi3)

        # Minus strand version of above
        # [0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2]
        window_size = 6
        win_start, win_end = g2.getTssCentredCoords(window_size)
        roi4 = ROI(g2, win_start, win_end)
        ROIs.append(roi4)

        # Test alternative chromosomes of interest
        # [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        window_size = 6
        win_start, win_end = g3.getTssCentredCoords(window_size)
        roi5 = ROI(g3, win_start, win_end)
        ROIs.append(roi5)

        # NOTE: No alternative chroms in test BED file to test for...

        return ROIs

    def test_count_BAM(self):
        """ make sure we get the expected counts from BAM reading """
        ROIs = self.setUpROIsForFiles()
        ROIs, rr = read_BAM('data/brca2-11_sorted.bam', ROIs)
        self.assertEqual(ROIs[0].counts, [0, 0, 0, 0, 0, 2, 2, 2, 2, 2])
        self.assertEqual(ROIs[1].counts, [0, 0, 0, 0, 0, 0, 2, 2, 2, 2])
        self.assertEqual(ROIs[2].counts, [0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2])
        self.assertEqual(ROIs[3].counts, [0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2])
        self.assertEqual(ROIs[4].counts, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

    def test_count_BED(self):
        """ make sure we get the expected counts from BED reading
        read_BED will often mess up the order of the ROIs so we'll need
        a clever method to index them to check them correctly """
        ROIs = self.setUpROIsForFiles()
        ROIs, rr = read_BED('data/brca2-11.bed', ROIs)

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
                self.assertEqual(roi.counts, 'Unknown')
