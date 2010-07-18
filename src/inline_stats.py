from __future__ import division

import numpy, cPickle
from cogent.maths.stats.util import NumberFreqs, Numbers

def dict_to_array(data):
    qual_char = dict([(ord(c)-66, c) for c in data])
    counts = numpy.zeros([len(data), len(data[c])], int)
    positions = numpy.arange(counts.shape[1])
    for i in range(counts.shape[0]):
        c = qual_char[i]
        counts[i].put(positions, data[c])
    return counts

class RunningStats(object):
    def __init__(self, data=None, in_file=None):
        super(RunningStats, self).__init__()
        assert data is not None or in_file is not None
        self._mean = None
        self._sd = None
        
        if in_file is not None:
            self.loadStats(in_file)
        else:
            self._counts = dict_to_array(data)
        
        self.length = self._counts.shape[1]
    
    def _get_mean(self):
        if self._mean is not None:
            return self._mean
        
        column_sums = self._counts.sum(axis=0)
        row_sums = self._counts.sum(axis=1)
        result = []
        for p, num in enumerate(column_sums):
            if num == 0:
                continue
            col_total = 0
            for qual, count in enumerate(self._counts[:, p]):
                col_total += ((2+qual) * count)
            result.append(col_total / num)
        self._mean = numpy.array(result)
        return self._mean
    
    Mean = property(_get_mean)

    def _get_sd(self):
        if self._sd is not None:
            return self._sd
        
        pos_means = self.Mean
        pos_sums = self._counts.sum(axis=0)
        result = []
        for pos in range(pos_means.shape[0]):
            mean = pos_means[pos]
            var = 0.0
            for qual in range(38):
                var += ((qual + 2 - mean)**2  * self._counts[qual, pos])
            result.append(numpy.sqrt(var / pos_sums[pos]))
        self._sd = numpy.array(result)
        return self._sd
    
    SD = property(_get_sd)
    
    def quantiles(self, quantiles_):
        """returns position-wise quantiles
        
        Uses Cogent's NumberFreqs and Numbers objects"""
        # this is very inefficient
        quals = [qual+2 for qual in range(38)]
        results = []
        for position in range(self._counts.shape[1]):
            column = NumberFreqs(data=dict(zip(quals,
                        self._counts[:, position])))
            column = Numbers(column.expand())
            results += [[column.quantile(q) for q in quantiles_]]
        
        return numpy.array(results).T
    
    def storeStats(self, out_file=None):
        """Store stats in pickle file for easy retrieval"""

        if out_file is None:
            out_file = 'stats.pkl'

        pklFile = open(out_file, 'wb')
        cPickle.dump(self._counts, pklFile)
        pklFile.close()

    def loadStats(self, in_file):
        """Initialise the object from previously saved data"""
        
        pklFile = file(in_file)
        self._counts = cPickle.load(pklFile)
        pklFile.close()

