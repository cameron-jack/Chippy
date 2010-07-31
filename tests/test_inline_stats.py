import numpy
import sys
sys.path.append('../src')

from cogent.util.unit_test import TestCase, main
from cogent.maths.stats.util import NumberFreqs, Numbers
from cogent.util.misc import remove_files

from light_seq import LightSeq
from inline_stats import RunningStats, quantile
from parse_fastq import FastqParser

def _to_freqs(values):
    """convenience func to produce a value, counts structure for checking
    quantile code"""
    data = {}
    for value in values:
        if value in data:
            data[value] += 1
        else:
            data[value] = 1
    return zip(*sorted(data.items()))

def _get_data_dict(length):
    """returns a dict for incrementing counts"""
    data = dict([(chr(i+66), [0]*length) for i in range(39)])
    return data

def RecordQuality(data):
    """records qutlity scores in data dict"""
    def call(qual):
        for i in range(len(qual)):
            data[qual[i]][i] += 1
    return call

class RunningStatsTest(TestCase):


    def test_save_load(self):
        """Making sure that saving / loading Running Stats data works"""
        seqs = [a for a in FastqParser('data/test.txt')]

        data = _get_data_dict(35)
        qualities = RecordQuality(data)
        maxSeqLength = 0
        for seq in seqs:
            qualities(seq.Quality)
            seqLength = len(seq)
            if seqLength > maxSeqLength:
                maxSeqLength = seqLength

        stats = RunningStats(data)
        # Also try saving / initialising the object from stored file
        pickle_file = 'data/test_stats.pkl'
        stats.storeStats(pickle_file)
        loaded_stats = RunningStats(in_file=pickle_file)
        self.assertFloatEqual(stats.Mean, loaded_stats.Mean)
        remove_files([pickle_file])

    def est_stats(self):
        """Test whether the mean and the std. deviation computed by RunningStats
        approximates real mean and std. deviation. Mean should be exact, however
        standard deviation is an approximation"""

        seqs = [a for a in FastqParser('data/test.txt')]
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

    def test_quantiles(self):
        """quantile funcs should match cogent"""
        seqs = [a for a in FastqParser('data/test.txt')]
        qualList = []
        data = _get_data_dict(35)
        qualities = RecordQuality(data)
        for seq in seqs:
            qualList.append(seq.getNumericQuality())
            qualities(seq.Quality)

        qualArray = numpy.array(qualList)
        for position in range(qualArray.shape[0]):
            freqs = Numbers(qualArray[0, :])
            for quant in 0.05, 0.5, 0.95:
                values, counts = _to_freqs(qualArray[0, :])
                got = quantile(numpy.array(values), numpy.array(counts), quant)
                expect = freqs.quantile(quant)
                self.assertFloatEqual(got, expect)


if __name__ == "__main__":
    main()