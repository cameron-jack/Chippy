from __future__ import division

import numpy, cPickle

_data = numpy.zeros((38, 75), int)

class RunningStats(object):
    """computes running mean and SD"""
    def __init__(self, length=1, in_file=None):
        super(RunningStats, self).__init__()
        self.length = length
        self.sumx = numpy.zeros(length, numpy.uint64)
        self.sumx2 = numpy.zeros(length, numpy.uint64)
        self.counts = numpy.zeros(length, numpy.uint32)

        if in_file is not None:
            self.loadStats(in_file)

    def __call__(self, vals):
        """update the sum, and sum of squared and counts when called"""
        for i, v in enumerate(vals):
            qual_index = v-2
            _data[qual_index, i] += 1
        
        # grow the length if needed
        valLength = len(vals)
        if valLength > self.length:
            self.sumx.resize(valLength)
            self.sumx2.resize(valLength)
            self.counts.resize(valLength)
            self.length = valLength

        index = 0
        for val in vals:
            self.counts[index] += 1
            self.sumx[index] += val
            self.sumx2[index] += (val*val)
            index += 1


    def _get_mean(self):
        """Compute the mean for each index"""
        mean = numpy.zeros(self.length, float)
        index = 0
        for val in self.sumx:
            mean[index] = (float(val) / self.counts[index])
            index += 1

        return mean
    
    def _get_array_based_mean(self):
        column_sums = _data.sum(axis=0)
        row_sums = _data.sum(axis=1)
        result = []
        for p, num in enumerate(column_sums):
            if num == 0:
                continue
            col_total = 0
            for qual, count in enumerate(_data[:, p]):
                col_total += ((2+qual) * count)
            result.append(col_total / num)
        return numpy.array(result)
    
    Mean = property(_get_mean)

    def _get_sd(self):
        """Compute the standard deviation for each index"""
        sd = numpy.zeros(self.length, float)
        index = 0
        for x2sum in self.sumx2:
            xsum = self.sumx[index]
            n = self.counts[index]
            sd[index] = numpy.sqrt((x2sum - (xsum*xsum / n))/n)
            index += 1

        return numpy.around(sd, 2)
    
    def _get_array_based_sd(self):
        col_means = self._get_array_based_mean()
        column_sums = _data.sum(axis=0)
        row_sums = _data.sum(axis=1)
        result = []
        for pos in range(col_means.shape[0]):
            mean = col_means[pos]
            var = 0.0
            for qual in range(38):
                var += ((mean - qual)**2  * _data[qual, pos])
            result.append(numpy.sqrt(var / column_sums[pos]))
        return numpy.array(result)
    
    SD = property(_get_sd)

    def storeStats(self, out_file=None):
        """Store stats in pickle file for easy retrieval"""

        if out_file is None:
            out_file = 'stats.pkl'

        pklFile = open(out_file, 'wb')
        cPickle.dump(self.sumx, pklFile)
        cPickle.dump(self.sumx2, pklFile)
        cPickle.dump(self.counts, pklFile)
        cPickle.dump(self.length, pklFile)
        pklFile.close()

    def loadStats(self, in_file):
        """Initialise the object from previously saved data"""

        pklFile = file(in_file)
        self.sumx = cPickle.load(pklFile)
        self.sumx2 = cPickle.load(pklFile)
        self.counts = cPickle.load(pklFile)
        self.length = cPickle.load(pklFile)
        pklFile.close()



