import numpy
import sys
sys.path.append('../src')

from cogent.util.unit_test import TestCase, main
from light_seq import LightSeq
from inline_stats import RunningStats
from parse_fastq import FastqParser



class RunningStatsTest(TestCase):


    def test_load(self):
        """Making sure that the loading of the Running Stats data works"""
        seqs = [a for a in FastqParser('data/test.txt', numeric_qual=False, make_seq=LightSeq)]
        stats = RunningStats()
        maxSeqLength = 0
        for seq in seqs:
            stats(seq.getNumericQuality())
            seqLength = len(seq.Seq)
            if seqLength > maxSeqLength:
                maxSeqLength = seqLength


        self.assertEqual(stats.counts[0], len(seqs))
        self.assertEqual(stats.length, maxSeqLength)

    def test_stats(self):
        """Test whether the mean and the std. deviation computed by RunningStats
        approximates real mean and std. deviation. Mean should be exact, however
        standard deviation is an approximation"""

        seqs = [a for a in FastqParser('data/test.txt', numeric_qual=False, make_seq=LightSeq)]
        stats = RunningStats()

        qualList = []
        for seq in seqs:
            stats(seq.getNumericQuality())
            qualList.append(seq.getNumericQuality())

        # compute actual values
        qualArray = numpy.array(qualList)
        mean = qualArray.mean(axis=0)
        sd = numpy.around(qualArray.std(axis=0), 2)

        self.assertEqual(stats.Mean, mean)
        self.assertFloatEqualRel(stats.SD, sd, eps=0.05)


if __name__ == "__main__":
    main()