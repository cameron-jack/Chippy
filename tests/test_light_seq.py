#!/usr/bin/env python

import sys
sys.path.append('../src')

from cogent.util.unit_test import TestCase, main
from light_seq import LightSeq

test_seq_name = 'GAPC_0015:6:1:1259:10413#0/1'
test_seq = 'AACACCAAACTTCTCCACCACGTGAGCTACAAAAG'
test_seq_qual = r'````Y^T]`]c^cabcacc`^Lb^ccYT\T\Y\WF'
test_seq_fasta = '>GAPC_0015:6:1:1259:10413#0/1\nAACACCAAACTTCTCCACCACGTGAGCTACAAAAG'
test_seq_fastq = '@GAPC_0015:6:1:1259:10413#0/1\nAACACCAAACTTCTCCACCACGTGAGCTACAAAAG\n+GAPC_0015:6:1:1259:10413#0/1\n````Y^T]`]c^cabcacc`^Lb^ccYT\T\Y\WF'

class LightSeqTest(TestCase):
    def test_load(self):
        """test the loading of the data into the sequence object"""

        sequence = LightSeq(Seq=test_seq, Name=test_seq_name, Quality=test_seq_qual)
        self.assertEqual(sequence.Seq, test_seq)
        self.assertEqual(sequence.Name, test_seq_name)
        self.assertEqual(sequence.Quality, test_seq_qual)

        sequence(test_seq, test_seq_name, test_seq_qual)
        self.assertEqual(sequence.Seq, test_seq)
        self.assertEqual(sequence.Name, test_seq_name)
        self.assertEqual(sequence.Quality, test_seq_qual)


    def test_numeric_conversion(self):
        """ test the numeric conversion of the quality scores in the sequence
        object"""

        sequence = LightSeq(Seq=test_seq, Name=test_seq_name, Quality=test_seq_qual)

        # Illumina schema allows for quality scores between 2 and 40 where 2 is
        # is considered a Read Segment Quality Control Indicator (towards the end
        # of the read
        qualScores = sequence.getNumericQuality(scheme='illumina')
        for score in qualScores:
            self.assertTrue(2 <= score <= 40)

        sequence(test_seq, test_seq_name, test_seq_qual)
        qualScores = sequence.getNumericQuality(scheme='illumina')
        for score in qualScores:
            self.assertTrue(2 <= score <= 40)



    def test_trim_seqs(self):
        """ test that the trimming of sequences is done atomically i.e. if a
        sequence is trimmed, so is the quality score."""

        sequence = LightSeq(Seq=test_seq, Name=test_seq_name, Quality=test_seq_qual)
        new_sequence = sequence[:-10]
        self.assertEqual(new_sequence.Seq, test_seq[:-10])
        self.assertEqual(new_sequence.Quality, test_seq_qual[:-10])

        sequence(test_seq, test_seq_name, test_seq_qual)
        new_sequence = sequence[:-10]
        self.assertEqual(new_sequence.Seq, test_seq[:-10])
        self.assertEqual(new_sequence.Quality, test_seq_qual[:-10])


    def test_toFasta(self):
        """ test fasta conversion"""

        sequence = LightSeq(Seq=test_seq, Name=test_seq_name, Quality=test_seq_qual)
        self.assertEqual(sequence.toFasta(), test_seq_fasta)

        sequence(test_seq, test_seq_name, test_seq_qual)
        self.assertEqual(sequence.toFasta(), test_seq_fasta)


    def test_toFastq(self):
        """ test fastq conversion"""

        # Quality score specified at creation of object
        sequence = LightSeq(Seq=test_seq, Name=test_seq_name, Quality=test_seq_qual)
        self.assertEqual(sequence.toFastq(), test_seq_fastq)

        # Quality score specified when creating fastq string
        sequence = LightSeq(Seq=test_seq, Name=test_seq_name)
        self.assertEqual(sequence.toFastq(test_seq_qual), test_seq_fastq)

        sequence(test_seq, test_seq_name, test_seq_qual)
        self.assertEqual(sequence.toFastq(test_seq_qual), test_seq_fastq)


if __name__ == "__main__":
    main()