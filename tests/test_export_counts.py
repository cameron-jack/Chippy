import sys
sys.path.append('../src')

import numpy, os

from cogent.util.unit_test import TestCase, main
from cogent import LoadTable
from cogent.util.misc import remove_files
from export_counts import get_counts
from util import data_dir


stable_id_file = './data/stable_ids.txt'
ctl_path = './data'
ctl_file = './data/s_8-window_10000-chrY.npy'
ctl_lane = 8
trt_path = './data'
trt_file = './data/s_2-window_10000-chrY.npy'
trt_lane = 2
score_path = './data'
score_file = './data/mouse_aln_chrY-10000-mapscore-chrY.npy'
gene_indices = LoadTable(os.path.join(data_dir,
                    'mouse_gene_coords.txt'), sep='\t')

class ExportCountsTests(TestCase):
    def test_get_counts(self):
        """Just checking that the function returns without failing"""
        table = get_counts(stable_id_file, ctl_path, ctl_lane, trt_path,
                           trt_lane, score_path)

    def test_valid_counts(self):
        """Validate that writing a subset table to file, and then rereading it
        reproduces the original counts values"""
        expected_counts_ctl = numpy.load(ctl_file)
        expected_counts_trt = numpy.load(trt_file)
        expected_scores = numpy.load(score_file)
        table = get_counts(stable_id_file, ctl_path, ctl_lane, trt_path,
                           trt_lane, score_path)
        table.writeToFile('test_counts.txt', sep=',')
        read_table = LoadTable('test_counts.txt', sep=',')
        for row in read_table.getRawData():
            stable_id = row.pop(0)
            index = gene_indices.filtered(lambda x: x == stable_id,
                                          columns=['StableId']).getRawData(columns=['Index'])[0]
            mode = row.pop(0) # control or treatment?
            if mode is 'c':
                self.assertEqual(numpy.array(row), expected_counts_ctl[index])
            elif mode is 't':
                self.assertEqual(numpy.array(row), expected_counts_trt[index])
            elif mode is 'm':
                self.assertFloatEqual(numpy.array(row), expected_scores[index])

        remove_files(['test_counts.txt'])


if __name__ == "__main__":
    main()