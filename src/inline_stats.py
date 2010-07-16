from __future__ import division

import numpy, cPickle
from cogent.maths.stats.util import NumberFreqs, Numbers

def resized(array, length):
    # required because a.resize broken in latest numpy release
    new = numpy.zeros((array.shape[0], length), numpy.uint64)
    indices = range(array.shape[1])
    for i in range(array.shape[0]):
        new[i].put(indices, array[i])
    return new

class RunningStats(object):
    """computes running mean and SD"""
    def __init__(self, length=1, in_file=None):
        super(RunningStats, self).__init__()
        self.length = length
        self._mean = None
        self._sd = None
        self._counts = numpy.zeros((38, length), numpy.uint64)
        
        if in_file is not None:
            self.loadStats(in_file)

    def __call__(self, vals):
        """record counts of quality scores when called"""
        self._mean = None
        self._sd = None
        
        if len(vals) > self.length:
            self.length = len(vals)
            self._counts = resized(self._counts, self.length)
        
        for i, v in enumerate(vals):
            qual_index = v-2
            self._counts[qual_index, i] += 1
        
    
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
    
    def quantile(self, q):
        """returns position-wise quantiles
        
        Uses Cogent's NumberFreqs and Numbers objects"""
        quals = [qual+2 for qual in range(38)]
        results = []
        for position in range(self._counts.shape[1]):
            column = NumberFreqs(data=dict(zip(quals,
                        self._counts[:, position])))
            column = Numbers(column.expand())
            results += [column.quantile(q)]
        return results
    
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



