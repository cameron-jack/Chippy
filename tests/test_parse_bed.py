import sys
sys.path.extend(['..', '../src'])

from cogent.util.unit_test import TestCase, main

from chippy.parse.bed import MinimalBedParser, BedRep

class ParseBed(TestCase):
    def test_BedRep_brca2_11(self):
        """ ests the BedRep class - relies on parser being correct """
        bed_data = BedRep('data/brca2-11.bed')
        chrom_5 = bed_data.get_chrom_as_nparray('5')

        parser = MinimalBedParser('data/brca2-11.bed')
        line_num = 0
        for expected_line in parser:
            self.assertEqual(expected_line[1], chrom_5[line_num,0])
            self.assertEqual(expected_line[2] - expected_line[1], chrom_5[line_num,1])
            self.assertEqual(expected_line[5], chrom_5[line_num,2])
            line_num += 1

        # second lines are the same


    def test_parser_brca2_11(self):
        """ tests the BED parser """

        #chr_5   151341223       151341298       ENSMUSG00000041147_151341223_plus_0     37      +
        #chr_5   151341223       151341298       ENSMUSG00000041147_151341223_minus_0    37      -
        expected_start = 151341223
        expected_length = 75
        expected_name = 'ENSMUSG00000041147_151341223'
        expected_score = 37
        expected_strand = 1

        running_length = 0

        parser = MinimalBedParser('data/brca2-11.bed')
        for line in parser:
            track_chrom = line[0]
            track_start = line[1]
            track_end = line[2]
            track_name = line[3]
            track_score = line[4]
            track_strand = line[5]

            # compare
            self.assertEqual('5', track_chrom)
            self.assertEqual(expected_start, track_start)
            self.assertEqual(expected_start + expected_length, track_end)
            self.assertEqual(expected_score, track_score)
            self.assertEqual(expected_strand, track_strand)
            if expected_strand == 1:
                test_name = '%s_plus_%d' % (expected_name, running_length)
                self.assertEqual(test_name, track_name)
            else:
                test_name = '%s_minus_%d' % (expected_name, running_length)
                self.assertEqual(test_name, track_name)

            # once we've seen the negative strand, advance our position
            if expected_strand == -1:
                running_length += 75
                expected_start += 75
            # swap strandedness
            expected_strand = [1,-1][expected_strand == 1]

if __name__ == '__main__':
    main()

