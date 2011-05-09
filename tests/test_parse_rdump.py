import sys
sys.path.append('..')

import numpy
from cogent import LoadTable
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files

from chippy.parse.r_dump import remove_probesets, SimpleRdumpToTable, convert

_sample_dump = LoadTable(header=['ENSEMBL', 'probeset', 'exp'],
        rows=[['id1',"0|1|2","13.6|13.4|13.6"],
        ['id2',"3|1","9.9|13.6"], # this gene should be lost when filtered
        ['id3',"4|5","12.7|13.4"],
        ['id4',"6","13.4"],
        ['id5',"7|8|3","6.0|6.0|4.5"],
        ['id6',"9|10|11|12","5.4|6.8|6.6|6.2"],
        ['id8',"13","12.7"],
        ['id9',"14","12.7"],
        ['id10',"15","12.7"]])

class ExcludingProbesets(TestCase):
    """test that excluding probesets works correctly"""
    def test_convert(self):
        """correctly cast | separated strings"""
        to_int = convert(to_float=False) # casts to int
        self.assertEqual(to_int('1|2|3'), (1,2,3))
        # but when not strict, returns strings if cannot cast
        self.assertEqual(to_int('a1|b2|c3'), ('a1','b2','c3'))
        to_int = convert(to_float=False, strict=True) # still casts to int
        self.assertEqual(to_int('1|2|3'), (1,2,3))
        # but when strict, and cannot cast raises ValueError
        self.assertRaises(ValueError, to_int, 'a1|b2|c3')
    
    def test_no_probes_remaining(self):
        """remove_probesets returns None if all probesets removed"""
        row = ['someid', (123, 345), (10.0,20.0)]
        remove = (123, 345)
        new_row = remove_probesets(row, remove, 1, 2)
        self.assertEqual(new_row, None)
    
    def test_allowed_probes_remain(self):
        """remove_probesets returns only probesets mapping to one gene"""
        row = ['someid', (123, 345), (10.0,20.0)]
        remove = (345,)
        new_row = remove_probesets(row, remove, 1, 2)
        self.assertEqual(new_row, ['someid', (123,), (10.0,)])
        remove = (123,)
        new_row = remove_probesets(row, remove, 1, 2)
        self.assertEqual(new_row, ['someid', (345,), (20.0,)])
        row = ['someid', (123, 345, 678), (10.0, 20.0, 30.0)]
        remove = (123, 678)
        new_row = remove_probesets(row, remove, 1, 2)
        self.assertEqual(new_row, ['someid', (345,), (20.0,)])
    
    def test_probeset_clean_fails_invalid(self):
        """remove_probesets should fail if an invalid probeset is provided"""
        row = ['someid', (123, 345, 678), (10.0, 20.0, 30.0)]
        self.assertRaises(ValueError, remove_probesets, row, (879,), 1, 2)
    
    def test_parsing_rdump(self):
        """parsing an rdump with no filtering returns same number of rows"""
        tmp_file = 'sample_rdump.txt'
        _sample_dump.writeToFile(tmp_file, sep='\t')
        rdumped, rr = SimpleRdumpToTable(tmp_file, stable_id_label='ENSEMBL',
                probeset_label='probeset', exp_label='exp',
                allow_probeset_many_gene=True)
        self.assertEqual(rdumped.Shape[0], _sample_dump.Shape[0])
        self.assertEqual(type(rdumped[0, 'probeset'][0]), int)
        self.assertEqual(type(rdumped[0, 'exp'][0]), float)
        remove_files([tmp_file], error_on_missing=False)
    
    def test_filtered_rdump(self):
        """filtered rdumps should exclude probesets that map to many genes"""
        tmp_file = 'sample_rdump.txt'
        _sample_dump.writeToFile(tmp_file, sep='\t')
        rdumped, rr = SimpleRdumpToTable(tmp_file, stable_id_label='ENSEMBL',
                probeset_label='probeset', exp_label='exp',
                allow_probeset_many_gene=False)
        self.assertTrue(rdumped.Shape[0] < _sample_dump.Shape[0])
        redundant_probesets = set([1,3])
        good_probesets = set(range(16)) - redundant_probesets
        probesets = [p for grp in rdumped.getRawData('probeset') for p in grp]
        probesets = set(probesets)
        self.assertEqual(probesets & redundant_probesets, set())
        self.assertEqual(probesets, good_probesets)
        remove_files([tmp_file], error_on_missing=False)
    
    def test_rdump_parse_validate(self):
        """validate arg to SimpleRdumpToTable rdump works"""
        tmp_file = 'sample_rdump.txt'
        sub = LoadTable(header=['ENSEMBL', 'probeset', 'exp'],
                rows=[['id1',"0|1|2","13.6|13.4|13.6"],
                ['id2',"3|1","9.9|13.6"],
                ['id3',"4|5","12.7|13.4"]])
        table = _sample_dump.appended(None, sub)
        table.writeToFile(tmp_file, sep='\t')
        # validate will cause an exception to be raised when True
        self.assertRaises(RuntimeError, SimpleRdumpToTable, tmp_file,
                stable_id_label='ENSEMBL', probeset_label='probeset',
                exp_label='exp', validate=True)
        # but not when False
        rdump, rr = SimpleRdumpToTable(tmp_file, stable_id_label='ENSEMBL',
            probeset_label='probeset', exp_label='exp', validate=False)
        
        # a different number of probeset, exp values should cause an exception
        sub = LoadTable(header=['ENSEMBL', 'probeset', 'exp'],
                rows=[['id1',"0|1|2","13.6|13.4|13.6"],
                ['id2',"3|1","9.9|13.6"],
                ['id3',"4|5","12.7"]])
        sub.writeToFile(tmp_file, sep='\t')
        self.assertRaises(RuntimeError, SimpleRdumpToTable, tmp_file,
                stable_id_label='ENSEMBL', probeset_label='probeset',
                exp_label='exp', validate=True)
        remove_files([tmp_file], error_on_missing=False)

if __name__ == '__main__':
    main()

