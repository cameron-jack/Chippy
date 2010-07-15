import numpy, cPickle

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



