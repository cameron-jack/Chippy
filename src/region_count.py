import time
from numpy import zeros, uint32, uint64, save, load

class RegionCounts(object):
    """records sequence read counts for a genomic region"""
    def __init__(self, length=None, filename=None):
        super(RegionCounts, self).__init__()
        
        assert length or filename, 'RegionCounts requires either'\
                                         ' a length or filename'
        if filename:
            self.counts = load(filename)
        else:
            self.counts = zeros((2, length), uint32)
    
    def addRead(self, start, end, strand):
        """add counts for range start, end"""
        
        self.counts[strand, start: end] += 1
    
    def getCounts(self, start, end, strand=None):
        """returns counts for the named span"""
        segment = self.counts[:, start: end]
        if strand is None:
            counts = segment.sum(axis=0)
        else:
            counts = segment[strand]
        return counts
    
    def save(self, filename):
        """saves as in numpy binary format"""
        save(filename, self.counts)
    

