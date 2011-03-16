import sys
sys.path.append('../src')

from cogent.util.unit_test import TestCase, main
from cogent import LoadTable
from cogent.parse.bowtie import BowtieOutputParser, BowtieToTable

fname = 'data/bowtie_output.txt'
data = [['GAPC_0015:6:1:1283:11957#0/1', '-', 'Mus', 66047927, 'TGTATATATAAACATATATGGAAACTGAATATATATACATTATGTATGTATATATGTATATGTTATATATACATA', 'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII', 0, ['55:A>G', '64:C>A']],
        ['GAPC_0015:6:1:1394:18813#0/1', '+', 'Mus', 77785518, 'ATGAAATTCCTAGCCAAATGGATGGACCTGGAGGGCATCATC', 'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII', 447, []],
        ['GAPC_0015:6:1:1560:18056#0/1', '+', 'Mus', 178806665, 'TAGATAAAGGCTCTGTTTTTCATCATTGAGAAATTGTTATTTTTCTGATGTTATA', 'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII', 0, ['9:T>G']],
        ['GAPC_0015:6:1:1565:19849#0/1', '+', 'Mus', 116516430, 'ACCATTTGCTTGGAAAATTGTTTTCCAGCCTTTCACTCTGAG', 'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII', 141, []],
        ['GAPC_0015:6:1:1591:17397#0/1', '-', 'Mus', 120440696, 'TCTAAATCTGTTCATTAATTAAGCCTGTTTCCATGTCCTTGGTCTTAAGACCAATCTGTTATGCGGGTGTGA', 'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII', 0, ['70:A>C', '71:G>T']]]

class BowtieOutputTest(TestCase):

    def test_parsing(self):
        """make sure that the bowtie output file is parsed properly"""

        parser = BowtieOutputParser(fname)
        header = parser.next()
        index = 0
        for row in parser:
            self.assertEqual(row, data[index])
            index += 1

    def test_psl_to_table(self):
        """make sure that the table is built without any hassle"""
        table = BowtieToTable(fname)

    def test_getting_seq_coords(self):
        """get correct information from the table"""
        table = BowtieToTable(fname)
        index = 0
        for row in table:
            query_name = row['Query Name']
            strand_direction = row['Strand Direction']
            query_offset = row['Offset']
            self.assertEqual(query_name, data[index][0])
            self.assertEqual(strand_direction, data[index][1])
            self.assertEqual(query_offset, data[index][3])
            index += 1

if __name__ == "__main__":
    main()