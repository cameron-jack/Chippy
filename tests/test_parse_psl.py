import sys
sys.path.append('../src')

from cogent.util.unit_test import TestCase, main
from cogent import LoadTable

from parse_psl import make_header, MinimalPslParser, PslToTable

fname = 'data/result_test_500_1.psl'

class Test(TestCase):

    def test_header(self):
        """should return correct header"""
        expect = ['match', 'mis-match', 'rep. match', "N's", 'Q gap count',
            'Q gap bases', 'T gap count', 'T gap bases', 'strand', 'Q name',
            'Q size', 'Q start', 'Q end', 'T name', 'T size', 'T start',
            'T end', 'block count', 'blockSizes', 'qStarts', 'tStarts']

        parser = MinimalPslParser(fname)
        version = parser.next()
        header = parser.next()
        self.assertEqual(header, expect)
        for row in parser:
            print row

    def test_psl_to_table(self):
        table = PslToTable(fname)

    def test_getting_seq_coords(self):
        """get correct sequence coordinates to produce a trimmed sequence"""
        table = PslToTable(fname)
        for row in table:
            query_name = row["Q name"]
            query_strand = row["strand"]
            q_start = row["Q start"]
            #print query_name, q_start, query_strand


if __name__ == "__main__":
    main()