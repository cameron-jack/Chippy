import sys
sys.path.extend(['../src', '..'])
import warnings
warnings.filterwarnings('ignore',
        "Not using MPI as mpi4py not found")

from chippy.core.region_of_interest import ROI
from cogent.util.unit_test import TestCase
from chippy.core.read_count import read_BAM, read_BED, read_BEDgraph
from chippy.express.db_schema import Gene

from chippy.util.definition import PLUS_STRAND, MINUS_STRAND, NULL_STRAND
import numpy

class MinimalRegionCountTests(TestCase):
    """
        Tests for: read_BAM, read_BED, read_BEDgraph
        Need test for read_wiggle
    """

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
        win_start, win_end = g1.getTssWindowCoords(window_radius, window_radius)
        roi1 = ROI(g1, win_start, win_end)
        ROIs.append(roi1)

        window_radius = 5
        # [0, 0, 0, 0, 0, 0, 2, 2, 2, 2] due to reads finishing 1 early
        win_start, win_end = g2.getTssWindowCoords(window_radius, window_radius)
        roi2 = ROI(g2, win_start, win_end)
        ROIs.append(roi2)

        # perfectly overlapping, and larger than roi1
        # This tests that the readers can handle overlapping ROIs
        # It also shows that the window shifts properly
        # [0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2]
        window_radius = 6
        win_start, win_end = g1.getTssWindowCoords(window_radius, window_radius)
        roi3 = ROI(g1, win_start, win_end)
        ROIs.append(roi3)

        # Minus strand version of above
        # [0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2]
        window_radius = 6
        win_start, win_end = g2.getTssWindowCoords(window_radius, window_radius)
        roi4 = ROI(g2, win_start, win_end)
        ROIs.append(roi4)

        # Test alternative chromosomes of interest
        # [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        window_radius = 6
        win_start, win_end = g3.getTssWindowCoords(window_radius, window_radius)
        roi5 = ROI(g3, win_start, win_end)
        ROIs.append(roi5)

        # NOTE: No alternative chroms in test BED file to test for...

        roi_ids = tuple([roi1.unique_id, roi2.unique_id,
                         roi3.unique_id, roi4.unique_id, roi5.unique_id])
        expected_arrays = tuple([
                numpy.array([0., 0., 0., 0., 0., 2., 2., 2., 2., 2.],
                        dtype=numpy.float32),
                numpy.array([0., 0., 0., 0., 0., 0., 2., 2., 2., 2.],
                        dtype=numpy.float32),
                numpy.array([0., 0., 0., 0., 0., 0., 2., 2., 2., 2., 2., 2.],
                        dtype=numpy.float32),
                numpy.array([0., 0., 0., 0., 0., 0., 0., 2., 2., 2., 2., 2.],
                        dtype=numpy.float32),
                numpy.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                        dtype=numpy.float32),
            ])

        return ROIs, roi_ids, expected_arrays

    def test_count_BAM(self):
        """
            Make sure we get the expected counts from BAM reading.
            Negative stand expectations are offset here because that's
            what's in the data
        """
        ROIs, roi_ids, expected = self.setUpROIsForFiles()
        ROIs, num_tags, num_bases, mapped_tags =\
                read_BAM('data/brca2-11_sorted.bam', ROIs)

        self.assertEqual(num_tags, 126) # total tags
        self.assertEqual(num_bases, 9450) # total bases
        self.assertEqual(mapped_tags, 126) # total experimental read tags

        for roi in ROIs:
            if roi.unique_id == roi_ids[0]:
                self.assertEqual(roi.counts, expected[0])
            elif roi.unique_id == roi_ids[1]:
                self.assertEqual(roi.counts, expected[1])
            elif roi.unique_id == roi_ids[2]:
                self.assertEqual(roi.counts, expected[2])
            elif roi.unique_id == roi_ids[3]:
                self.assertEqual(roi.counts, expected[3])
            else:
                self.assertEqual(roi.counts, expected[4])

    def test_count_BED(self):
        """
            Make sure we get the expected counts from BED reading
            read_BED will often mess up the order of the ROIs so we'll need
            a clever method to index them to check them correctly.
            Negative stand expectations are offset here because that's what's
            in the data.
        """
        ROIs, roi_ids, expected = self.setUpROIsForFiles()
        ROIs, num_tags, num_bases, mapped_tags = read_BED('data/brca2-11.bed', ROIs)

        self.assertEqual(num_tags, 126.0) # number of bed lines
        self.assertEqual(num_bases, 9450.0) # sum of all counts
        self.assertEqual(mapped_tags, 126) # total experimental read tags

        for roi in ROIs:
            if roi.unique_id == roi_ids[0]:
                self.assertEqual(roi.counts, expected[0])
            elif roi.unique_id == roi_ids[1]:
                self.assertEqual(roi.counts, expected[1])
            elif roi.unique_id == roi_ids[2]:
                self.assertEqual(roi.counts, expected[2])
            elif roi.unique_id == roi_ids[3]:
                self.assertEqual(roi.counts, expected[3])
            else:
                self.assertEqual(roi.counts, expected[4])

    def test_count_BEDgraph(self):
        """
            Make sure we get the expected counts from BEDgraph files.
        """
        ROIs, roi_ids, expected = self.setUpROIsForFiles()
        ROIs, num_tags, num_bases, mapped_tags = read_BEDgraph('data/brca2-11.bedgraph', ROIs)
        self.assertEqual(num_tags, 126) # tag count estimate
        self.assertEqual(num_bases, 9450) # sum of all counts
        self.assertEqual(mapped_tags, 126) # total experimental read tags

        for roi in ROIs:
            if roi.unique_id == roi_ids[0]:
                self.assertEqual(roi.counts, expected[0])
            elif roi.unique_id == roi_ids[1]:
                self.assertEqual(roi.counts, expected[1])
            elif roi.unique_id == roi_ids[2]:
                self.assertEqual(roi.counts, expected[2])
            elif roi.unique_id == roi_ids[3]:
                self.assertEqual(roi.counts, expected[3])
            else:
                self.assertEqual(roi.counts, expected[4])

