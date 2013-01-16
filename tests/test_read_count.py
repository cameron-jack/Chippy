import sys
sys.path.extend(['../src', '..'])
import warnings
warnings.filterwarnings('ignore',
        "Not using MPI as mpi4py not found")

from cogent.util.unit_test import TestCase, main
from chippy.core.count_tags import ROI
from chippy.core.read_count import read_BAM, read_BED

from chippy.util.definition import PLUS_STRAND, MINUS_STRAND, NULL_STRAND

class MinimalRegionCountTests(TestCase):
    """
    tests for: read_BAM, read_BED,ROI.
    """

    def setUpROIs(self):
        """
        In brca2-11.sam if we make TSS= 151341224 and window size = 100.
        Then we should get 124->223 at 0 and 224->324 at +2.
        If we take TSS= 151345949 and window size = 100, then we should get
        849->949 at +2 and 950->1049 at 0, in 1-offset space.
        """
        window_size = 100
        ROIs = []

        TSS_1 = 151341224
        window_start = TSS_1 - window_size
        window_end = TSS_1 + window_size
        roi = ROI('chr_5', 'roi_1', 1, TSS_1, window_start, window_end,
                PLUS_STRAND)
        ROIs.append(roi)

        TSS_2 = 151345949
        window_start = TSS_1 - window_size
        window_end = TSS_1 + window_size
        roi = ROI('chr_5', 'roi_2', 2, TSS_2, window_start, window_end,
                PLUS_STRAND)
        ROIs.append(roi)

        return ROIs

    def test_count_BAM(self):
        """ make sure we get the expected counts from BAM reading """
        ROIs = self.setUpROIs()
        ROIs, rr = read_BAM('data/brca2-11_sorted.bam', ROIs)
        for roi in ROIs:
            print roi.counts

    def test_count_BED(self):
        """ make sure we get the expected counts from BAM reading """
        ROIs = self.setUpROIs()
        ROIs, rr = read_BED('data/brca2-11.bed', ROIs)
        for roi in ROIs:
            print roi.counts
