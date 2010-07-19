from __future__ import division

import numpy, cPickle

def value_at_expanded_index(values, counts, index):
    """return the value at index from expansion, using counts, of values"""
    num = 0
    cumsum = counts.cumsum()
    
    for i in range(cumsum.shape[0]):
        if cumsum[i] > index:
            return values[i]

def quantile(values, counts, quant):
    index = quant * (sum(counts)-1)
    lo = int(numpy.floor(index))
    hi = int(numpy.ceil(index))
    diff = index - lo
    lo_val = value_at_expanded_index(values, counts, lo)
    if diff != 0:
        hi_val = value_at_expanded_index(values, counts, hi)
    else:
        hi_val = 0
    stat = (1-diff) * lo_val + diff * hi_val
    return stat

def qual_counts_dict_to_array(data, scheme='illumina'):
    adjust = [33, 66][scheme == 'illumina']
    qual_char = dict([(ord(c)-adjust, c) for c in data])
    counts = numpy.zeros([len(data), len(data[c])], numpy.uint64)
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
            self._counts = qual_counts_dict_to_array(data)
        
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
            for qual in range(39):
                var += ((qual + 2 - mean)**2  * self._counts[qual, pos])
            result.append(numpy.sqrt(var / pos_sums[pos]))
        self._sd = numpy.array(result)
        return self._sd
    
    SD = property(_get_sd)
    
    def quantiles(self, quantiles_):
        """returns position-wise quantiles"""
        quals = numpy.array([qual+2 for qual in range(39)])
        results = []
        for position in range(self._counts.shape[1]):
            column = self._counts[:, position]
            indices = [i for i in range(len(column)) if column[i] != 0]
            qual = quals.take(indices)
            counts = column.take(indices)
            results += [[quantile(qual, counts, q) for q in quantiles_]]
        
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

