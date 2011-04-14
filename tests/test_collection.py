import sys
sys.path.append('..')

import numpy
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files

from chippy.core.collection import RegionCollection


class CollectionTests(TestCase):
    
    def check_write_load(self, orig_data):
        orig = RegionCollection(**orig_data)
        orig.writeToFile('test_data')
        recovered = RegionCollection(filename='test_data')
        for key in orig_data:
            self.assertEqual(getattr(recovered, key), orig_data[key])
        
        remove_files(['test_data'], error_on_missing=False)
    
    def test_save_load(self):
        """correctly save and reload a collection"""
        orig_data = dict(counts=['abcd', 'efg'],
                ranks=[0,1], labels=['A', 'B'], info={'args': 'some args'})
        self.check_write_load(orig_data)
    
    def test_consistency_check_on_create(self):
        """attributes need to have same number of elements"""
        orig_data = dict(counts=['abcd', 'efg'],
                ranks=[0,1,2], labels=['A', 'B'], info={'args': 'some args'})
        
        self.assertRaises(AssertionError, RegionCollection,
                        **orig_data)
        
        orig_data['ranks'] = [0,1]
        orig_data['labels'] = ['A']
        self.assertRaises(AssertionError, RegionCollection,
                        **orig_data)
    
    def test_grouping_data(self):
        """should correctly group data without a filter"""
        expect_two = dict(counts=[[0,1], [2,3], [4,5], [6,7], [8,9]],
            ranks=[[0,1], [2,3], [4,5], [6,7], [8,9]],
            labels=[['0','1'], ['2','3'], ['4','5'], ['6','7'], ['8','9']])
        
        coll = RegionCollection(**expect_two)
        
        grouped = coll.getGrouped(2)
        expect_counts_ranks = [[[0,1], [2,3]], [[4,5], [6,7]]]
        expect_labels = [[['0','1'], ['2','3']], [['4','5'], ['6','7']]]
        self.assertEqual(grouped[0].tolist(), expect_counts_ranks)
        self.assertEqual(grouped[1].tolist(), expect_counts_ranks)
        self.assertEqual(grouped[2].tolist(), expect_labels)
    
    def test_filtered_normalised(self):
        """filtered with normalised cutoff should work"""
        data = numpy.random.randint(0, 20, size=100)
        mat = data.reshape((10,10))
        # make two outlier sequences, the 3rd, the 5th
        mat[3,0] = 1000
        mat[5,9] = 500
        labels = map(str, range(10))
        coll = RegionCollection(counts=mat, labels=labels)
        new = coll.filteredNormalised(cutoff=3.0)
        expect_indices = [0,1,2,4,6,7,8,9]
        self.assertEqual(new.counts.shape[0], 8)
        for i, index in enumerate(expect_indices):
            self.assertEqual(new.counts[i], coll.counts[index])
            self.assertEqual(new.labels[i], coll.labels[index])
    
    def test_filtered_grouping_data(self):
        """should correctly filter records"""
        input_two = dict(counts=[[0,1], [2,3], [4,5], [6,7], [8,9]],
            ranks=[[0,1], [2,3], [4,5], [6,7], [8,9]],
            labels=[['0','1'], ['2','3'], ['4','5'], ['6','7'], ['8','9']])
        
        expect_counts_ranks = [[[0,1], [4,5]], [[6,7], [8,9]]]
        expect_labels = [[['0','1'], ['4','5']], [['6','7'], ['8','9']]]
        
        func = lambda x: x.min() != 2 and x.max() != 3
        
        coll = RegionCollection(**input_two)
        grouped = coll.getGrouped(2, filtered=func)
        self.assertEqual(grouped[0].tolist(), expect_counts_ranks)
        self.assertEqual(grouped[1].tolist(), expect_counts_ranks)
        self.assertEqual(grouped[2].tolist(), expect_labels)
    
    def test_only_counts(self):
        """correctly handle only being provided with counts data"""
        input_two = dict(counts=[[0,1], [2,3], [4,5], [6,7], [8,9]])
        
        # check write/read
        self.check_write_load(input_two)
        
        # check filter works
        expect_counts = [[[0,1], [4,5]], [[6,7], [8,9]]]
        expect_ranks_labels = None
        
        func = lambda x: x.min() != 2 and x.max() != 3
        
        coll = RegionCollection(**input_two)
        
        grouped = coll.getGrouped(2, filtered=func)
        self.assertEqual(grouped[0].tolist(), expect_counts)
        self.assertEqual(grouped[1].tolist(), expect_ranks_labels)
        self.assertEqual(grouped[2].tolist(), expect_ranks_labels)
    
    def test_itergroups(self):
        """should correctly generate groups of counts, ranks etc .."""
        input_two = dict(counts=[[0,1], [2,3], [4,5], [6,7], [8,9]],
            ranks=[[0,1], [2,3], [4,5], [6,7], [8,9]],
            labels=[['0','1'], ['2','3'], ['4','5'], ['6','7'], ['8','9']])
        
        coll = RegionCollection(**input_two)
        
        expected_labels = [[['0','1'], ['2','3']], [['4','5'], ['6','7']]]
        i = 0
        for c, r, l in coll.itergroups(group_size=2):
            
            self.assertEqual(c, r)
            self.assertEqual(l, expected_labels[i])
            i += 1
        
        expect_two = input_two.copy()
        expect_two['counts'].remove([2,3])
        expect_two['ranks'].remove([2,3])
        expect_two['labels'].remove(['2','3'])
        
        expected_labels=[[['0','1'], ['4','5']], [['6','7'], ['8','9']]]
        func = lambda x: x.min() != 2 and x.max() != 3
        i = 0
        for c, r, l in coll.itergroups(group_size=2, filtered=func):
            self.assertEqual(c, r)
            self.assertEqual(l, expected_labels[i])
            i += 1
    
    def test_filtered_by_label(self):
        """return correct subset by label"""
        input = dict(counts=[[0,1], [2,3], [4,5], [6,7], [8,9]],
            ranks=[0, 1, 2, 3, 4],
            labels=['0','1', '2','3', '4'])
        coll = RegionCollection(**input)
        new = coll.filteredByLabel('1') # one label
        self.assertEqual(new.labels, ['1'])
        self.assertEqual(new.counts.tolist(), [[2,3]])
        self.assertEqual(new.ranks.tolist(), [1])
        
        new = coll.filteredByLabel(['1']) # one label
        self.assertEqual(new.labels, ['1'])
        self.assertEqual(new.counts.tolist(), [[2,3]])
        self.assertEqual(new.ranks.tolist(), [1])
        
        new = coll.filteredByLabel(['0', '2', '4']) # 3 disjoint labels
        self.assertEqual(new.labels, ['0', '2', '4'])
        self.assertEqual(new.counts.tolist(), [[0,1], [4,5], [8,9]])
        self.assertEqual(new.ranks.tolist(), [0, 2, 4])
        
        # without ranks
        input.pop('ranks')
        coll = RegionCollection(**input)
        new = coll.filteredByLabel('1') # one label
        self.assertEqual(new.labels, ['1'])
        self.assertEqual(new.counts.tolist(), [[2,3]])
        self.assertEqual(new.ranks, None)
    
    def test_filtered_by_label_fails(self):
        """filtered by label faisl if no labels"""
        input = dict(counts=[[0,1], [2,3], [4,5], [6,7], [8,9]])
        coll = RegionCollection(**input)
        self.assertRaises(RuntimeError, coll.filteredByLabel, '1')
    

if __name__ == '__main__':
    main()
