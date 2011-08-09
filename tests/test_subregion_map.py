import sys
sys.path.append('../src')

import numpy

from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files
from subregion_map import SubregionMap

start_sites_Y = [133852, 155091, 234230, 346985, 582202, 623020, 796225,
                 1426357, 1919568, 2086608, 2118067, 2160042]
coords = [(tss-10000, tss+10000) for tss in start_sites_Y]
Y_length = 15902555
counts_file = 'data/chrY.npy'
scores_file = 'data/chrY_score.npy'


class SubRegionMapTest(TestCase):

    def test_init(self):
        """Test the appropriate initialisation of the object"""

        mapscore = SubregionMap(counts_file, 'Y', Y_length, start_sites_Y, 75)
        self.assertEqual(mapscore.chrom, 'Y')
        self.assertEqual(mapscore.coords, coords)
        self.assertEqual(mapscore.window_size, 10000)

    def test_getMapScores(self):
        """ Test the correct computation of the mappability scores"""

        mapscore = SubregionMap(counts_file, 'Y', Y_length, start_sites_Y, 75)
        expected = numpy.load(scores_file)
        self.assertFloatEqual(mapscore.getMappabilityScore(), expected)

    def test_saveMapScores(self):
        """ Test that the scores are saved appropriately"""

        file_name = 'mappability_score'
        mapscore = SubregionMap(counts_file, 'Y', Y_length, start_sites_Y, 75)
        mapscore.saveMappabilityScore(file_name)
        got = numpy.load(file_name+'.npy')
        expect = numpy.load(scores_file)
        self.assertFloatEqual(got, expect)
        remove_files([file_name+'.npy'])

if __name__ == "__main__":
    main()
