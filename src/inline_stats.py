import numpy

class RunningStats(object):
    """computes running mean and SD"""
    def __init__(self, length):
        super(RunningStats, self).__init__()
        self.length = length
        self.mean = numpy.zeros(length, float)
        self.meansqr = numpy.zeros(length, float)
        self.counts = numpy.zeros(length, float)
    
    def __call__(self, val):
        self.mean += val
        self.meansqr += (val*val)
        self.counts += 1
    
    def _get_mean(self):
        """docstring for _get_mean"""
        return self.mean / self.counts
    
    Mean = property(_get_mean)
    
    def _get_sd(self):
        return self.meansqr / self.counts
    
    SD = property(_get_sd)
    

if __name__ == "__main__":
    vals = range(10)
    stats = RunningStats(1)
    for val in vals:
        stats(val)
    
    print stats.Mean
    print stats.SD