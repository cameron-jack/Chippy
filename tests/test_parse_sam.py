import sys
sys.path.extend(['..', '../src'])

from cogent.util.unit_test import TestCase, main

from chippy.parse.sam import MinimalSamParser

class ParseSam(TestCase):
    def test_parser_brca2_11(self):
        """correctly parse strand, length location data for BRCA2 exon 11"""
        parser = MinimalSamParser('data/brca2-11.sam')
        mapped = dict([(line[0], line)
                        for line in parser if not type(line) == dict])
        gene_start = 151341223
        for read in mapped:
            if 'plus' in read:
                strand = 1
            else:
                strand = -1
            record = mapped[read]
            self.assertEqual(record[1], strand)
            # sam counting is 1-based, from the samtools manual
            # 4 POS 1-based leftmost POSition/coordinate of clipped sequence
            # HOWEVER, our default converter makes all numbers zero-based
            mapped_start = record[3]
            expect_start = int(read.split('_')[-1]) + gene_start
            self.assertEqual(expect_start, mapped_start)
            # all mapped reads are 75 long
            self.assertEqual(record[5], 75)
    

if __name__ == '__main__':
    main()
