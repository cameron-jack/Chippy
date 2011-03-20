import sys
sys.path.append('../src')

from cogent.parse.bowtie import BowtieOutputParser, BowtieToTable
from cogent import LoadSeqs, DNA

from cogent.util.unit_test import TestCase, main
from minimal_reads import mapped_coords, make_chrom_coord_table

class ReadingBowtieOutput(TestCase):
    def test_mapped_seqs_match_known(self):
        """bowties coordinates for brca2 should match those from ensembl"""
        parser = BowtieOutputParser('data/brca2-11.map', row_converter=None)
        header = parser.next()
        index = 0
        mapped = dict([(row[0], row) for row in parser])
        
        ref = LoadSeqs('data/brca2-11.fasta', moltype=DNA, aligned=False)
        fwd = ref.Seqs[0]
        
        gene_start = None
        for read in mapped:
            record = mapped[read]
            # should be able to correctly obtain the brca2 seq
            from_start = int(read.split('_')[-1])
            if gene_start is None:
                gene_start = int(read.split('_')[-3])
            mapped_location = int(record[3])
            self.assertTrue(gene_start <= mapped_location)
            
            length = len(record[4])
            start = int(record[3])
            rel_start = start - gene_start
            mapped_seq = record[4]
            expected_seq = fwd[rel_start: rel_start+length]
            self.assertEqual(str(mapped_seq), str(expected_seq))
            
    def test_read_table_correctly_stores_coords(self):
        """parsing the .map output for brca2 exon11 produces correct seq coords"""
        seq = LoadSeqs('data/brca2-11.fasta', moltype=DNA, aligned=False)
        seq = seq.Seqs[0]
        gene_start = 151341223
        coords = mapped_coords('data/brca2-11.map', 1000, False)
        table = make_chrom_coord_table(coords, 'chr5')
        self.assertEqual(table.getDistinctValues('length'), set([75]))
        self.assertEqual(table.getDistinctValues('freq'), set([1]))
        plus = sorted(table.filtered('strand == 1').getRawData('start'))
        minus = sorted(table.filtered('strand == -1').getRawData('start'))
        self.assertTrue(gene_start in plus)
        self.assertTrue(gene_start in minus)
        for i in range(1, len(plus)):
            self.assertEqual(plus[i]-plus[i-1], 75)
        
        for i in range(1, len(minus)):
            self.assertEqual(minus[i]-minus[i-1], 75)
        
        frags = []
        for start in plus:
            start = start - gene_start
            end = start + 75
            frags.append(str(seq)[start:end])
        # make sure sequences have same length
        self.assertTrue(''.join(frags) == str(seq)[:len(plus)*75])
        frags = []
        for start in minus:
            start = start - gene_start
            end = start + 75
            frags.append(str(seq)[start:end])
        
        self.assertTrue(''.join(frags) == str(seq)[:len(plus)*75])
        
    

if __name__ == '__main__':
    main()
