import sys
sys.path.append('..')

import numpy
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files

from chippy.core.collection import RegionCollection, column_sum, \
        tchebysheff_upper

class UtilityFuncTests(TestCase):
    def test_sum_column(self):
        """correctly sum the columns of a 2D array"""
        data = numpy.array([range(4), range(4,8)])
        summed = column_sum(data)
        expect = numpy.array([4, 6, 8, 10])
        self.assertEqual(summed, expect)
    
    def test_tchebysheff_upper(self):
        """correctly compute one-sided Tchebysehff inequality"""
        p = 0.05
        k = tchebysheff_upper(p)
        self.assertFloatEqual(1.0/(1+k**2), p)

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
    
    def test_filtered_tchebysheff(self):
        """filtered with Tchebysheff cutoff should work"""
        data = numpy.random.randint(0, 20, size=100)
        mat = data.reshape((10,10))
        # make two outlier sequences, the 3rd, the 5th
        mat[3,0] = 1000
        mat[5,9] = 500
        labels = map(str, range(10))
        coll = RegionCollection(counts=mat, labels=labels)
        new = coll.filteredTchebysheffUpper(p=0.1)
        expect_indices = [0,1,2,4,6,7,8,9]
        self.assertEqual(new.counts.shape[0], 8)
        for i, index in enumerate(expect_indices):
            self.assertEqual(new.counts[i], coll.counts[index])
            self.assertEqual(new.labels[i], coll.labels[index])
    
    def test_filtered_tchebysheff_raises(self):
        """raises an exception if invalid cutoff provided"""
        data = dict(counts=[[0,1], [2,3], [4,5], [6,7], [8,9]])
        coll = RegionCollection(**data)
        self.assertRaises(RuntimeError, coll.filteredTchebysheffUpper, -0.1)
        self.assertRaises(RuntimeError, coll.filteredTchebysheffUpper, 2.0)
    
    def test_filtered_data(self):
        """should correctly filter records"""
        input_two = dict(counts=[[0,1], [2,3], [4,5], [6,7], [8,9]],
            ranks=[0, 1, 2, 3, 4],
            labels=['0', '1', '2', '3', '4'])
        
        expect_counts = [[0,1], [4,5], [6,7], [8,9]]
        expect_ranks = [0, 2, 3, 4]
        expect_labels = ['0', '2', '3', '4']
        
        func = lambda x: x.min() != 2 and x.max() != 3
        
        coll = RegionCollection(**input_two)
        self.assertEqual(coll.N, 5)
        new = coll.filtered(func)
        self.assertEqual(new.N, 4)
        self.assertEqual(new.counts.astype(int).tolist(), expect_counts)
        self.assertEqual(new.ranks.astype(int).tolist(), expect_ranks)
        self.assertEqual(new.labels.tolist(), expect_labels)
    
    def test_only_counts(self):
        """correctly handle only being provided with counts data"""
        input_two = dict(counts=[[0,1], [2,3], [4,5], [6,7], [8,9]])
        
        # check write/read
        self.check_write_load(input_two)
        
        # check filter works
        expect_counts = [[[0,1], [4,5]], [[6,7], [8,9]]]
        
        func = lambda x: x.min() != 2 and x.max() != 3
        
        coll = RegionCollection(**input_two)
        coll = coll.filtered(func)
        got_counts, got_ranks, got_labels = coll.getGrouped(2)
        self.assertEqual(got_counts.tolist(), expect_counts)
        self.assertEqual(got_ranks, None)
        self.assertEqual(got_labels, None)
    
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
        """filtered by label fails if no labels"""
        input = dict(counts=[[0,1], [2,3], [4,5], [6,7], [8,9]])
        coll = RegionCollection(**input)
        self.assertRaises(RuntimeError, coll.filteredByLabel, '1')
    
    def test_asfloats(self):
        """conversion correctly returns instance with counts as floats"""
        input = dict(counts=[[0,1], [2,3], [4,5], [6,7], [8,9]])
        expect = numpy.array(input['counts'], dtype=float)
        coll = RegionCollection(**input)
        float_coll = coll.asfloats()
        self.assertEqual(float_coll.counts, coll.counts)
    
    def test_total_counts(self):
        """total counts should be correct"""
        input = dict(counts=[[0,1], [2,3], [4,5], [6,7], [8,9]])
        coll = RegionCollection(**input)
        self.assertEqual(coll.Total, 45)
    
    def test_asfreqs(self):
        """should correctly convert counts to freqs"""
        input = dict(counts=[[0,1], [2,3], [4,5], [6,7], [8,9]])
        coll = RegionCollection(**input)
        freqs = coll.asfreqs()
        self.assertEqual(freqs.Total, 1.0)
    
    def test_transformed(self):
        """correctly return transform counts"""
        input = dict(counts=[[0,1], [2,3], [4,5], [6,7], [8,9]])
        coll = RegionCollection(**input)
        c, r = coll.transformed()
        self.assertEqual(c, [4,5])
        freqs = coll.asfreqs()
        c, r = freqs.transformed(counts_func=column_sum)
        self.assertEqual(c, [20./45., 25./45])
    
    def test_take(self):
        """take returns subset with correct attributes"""
        input = dict(counts=[[0,1], [2,3], [4,5], [6,7], [8,9]],
                info={'abcd': 'efg'}, labels=['0', '1', '2', '3', '4'])
        coll = RegionCollection(**input)
        sub = coll.take([0,1,2])
        self.assertEqual(sub.info, input['info'])
        self.assertEqual(sub.counts.tolist(), [[0,1], [2,3], [4,5]])
        self.assertEqual(sub.labels.tolist(), ['0', '1', '2'])
    
    def test_str(self):
        """correct string representation"""
        input = dict(counts=[[0,1], [2,3], [4,5], [6,7], [8,9]],
                info={'abcd': 'efg'}, labels=['0', '1', '2', '3', '4'])
        coll = RegionCollection(**input)
        self.assertEqual(str(coll),
            'RegionCollection(num_records=5; has_ranks=False; has_labels=True)')

if __name__ == '__main__':
    main()
