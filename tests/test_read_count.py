import sys
sys.path.extend(['../src', '..'])
import warnings
warnings.filterwarnings('ignore',
        "Not using MPI as mpi4py not found")

from cogent.util.unit_test import TestCase, main
from chippy.core.count_tags import ROI
from chippy.core.read_count import read_BAM, read_BED, add_counts_to_ROI

from chippy.util.definition import PLUS_STRAND, MINUS_STRAND, NULL_STRAND

class MinimalRegionCountTests(TestCase):
    """
    tests for: add_counts_to_ROI, read_BAM, read_BED,ROI.
    """

    def setUpROIsForSynthTests(self):
        """
        Minimal regions that should be easy to test
        TSS is in 0-based space
        """
        window_size = 5
        ROIs = []

        tss = 5
        # So the slice is [0:10], which is 10 long
        window_start = tss - window_size
        window_end = tss + window_size
        roi = ROI('1', '1', 1, tss, window_start, window_end,
                PLUS_STRAND)
        ROIs.append(roi)

        tss = 5
        # technically this includes genome position 0 which doesn't exist
        window_start = tss - window_size
        window_end = tss + window_size
        roi = ROI('1', '2', 2, tss, window_start, window_end,
                MINUS_STRAND)
        ROIs.append(roi)

        return ROIs

    def test_add_counts_to_ROI(self):
        """ If this function works, the rest should easy """
        ROIs = self.setUpROIsForSynthTests()

        # First case, perfect overlap, slightly inside the ROI
        entry_start = 1 # start at genome position 1
        entry_end = 7 # end and include genome position 7
        add_counts_to_ROI(ROIs[0], entry_start, entry_end)
        self.assertEqual(ROIs[0].counts, [0, 1, 1, 1, 1, 1, 1, 1, 0, 0])
        add_counts_to_ROI(ROIs[1], entry_start, entry_end)
        self.assertEqual(ROIs[1].counts, [0, 0, 1, 1, 1, 1, 1, 1, 1, 0])

        # test both addition to previous ROI and positive over-the-edge counts
        entry_start = 7
        entry_end = 12
        add_counts_to_ROI(ROIs[0], entry_start, entry_end)
        self.assertEqual(ROIs[0].counts, [0, 1, 1, 1, 1, 1, 1, 2, 1, 1])
        add_counts_to_ROI(ROIs[1], entry_start, entry_end)
        self.assertEqual(ROIs[1].counts, [1, 1, 2, 1, 1, 1, 1, 1, 1, 0])

        # test addition and negative over-the-edge counts
        entry_start = -2
        entry_end = 2
        add_counts_to_ROI(ROIs[0], entry_start, entry_end)
        self.assertEqual(ROIs[0].counts, [1, 2, 2, 1, 1, 1, 1, 2, 1, 1])
        add_counts_to_ROI(ROIs[1], entry_start, entry_end)
        self.assertEqual(ROIs[1].counts, [1, 1, 2, 1, 1, 1, 1, 2, 2, 1])

    def setUpROIsForFiles(self):
        """
        In brca2-11.sam if we make TSS= 151341224 and window size = 5.
        Then we should get 119->223 at 0 and 224->229 at +2.
        If we take TSS= 151345949 and window size = 5, then we should get
        945->949 at +2 and 950->954 at 0, in 1-offset space.
        """
        window_size = 5
        ROIs = []

        TSS_1 = 151341224
        window_start = TSS_1 - window_size
        window_end = TSS_1 + window_size
        roi = ROI('chr_5', 'roi_1', 1, TSS_1, window_start, window_end,
                PLUS_STRAND)
        ROIs.append(roi)

        TSS_2 = 151345949
        window_start = TSS_2 - window_size
        window_end = TSS_2 + window_size
        roi = ROI('chr_5', 'roi_2', 2, TSS_2, window_start, window_end,
                PLUS_STRAND)
        ROIs.append(roi)

        # perfectly overlapping, and larger than roi_2
        # This tests that the readers can handle overlapping ROIs
        window_size = 6
        TSS_3 = 151345949
        window_start = TSS_2 - window_size
        window_end = TSS_2 + window_size
        roi = ROI('chr_5', 'roi_3', 3, TSS_3, window_start, window_end,
            PLUS_STRAND)
        ROIs.append(roi)

        return ROIs

    def test_count_BAM(self):
        """ make sure we get the expected counts from BAM reading """
        ROIs = self.setUpROIsForFiles()
        ROIs, rr = read_BAM('data/brca2-11_sorted.bam', ROIs)
        self.assertEqual(ROIs[0].counts, [0, 0, 0, 0, 0, 2, 2, 2, 2, 2])
        self.assertEqual(ROIs[1].counts, [2, 2, 2, 2, 2, 2, 2, 2, 2, 2])
        self.assertEqual(ROIs[2].counts, [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2])

    def test_count_BED(self):
        """ make sure we get the expected counts from BAM reading """
        ROIs = self.setUpROIsForFiles()
        ROIs, rr = read_BED('data/brca2-11.bed', ROIs)
        self.assertEqual(ROIs[0].counts, [0, 0, 0, 0, 0, 2, 2, 2, 2, 2])
        self.assertEqual(ROIs[1].counts, [2, 2, 2, 2, 2, 2, 2, 2, 2, 2])
        self.assertEqual(ROIs[2].counts, [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2])