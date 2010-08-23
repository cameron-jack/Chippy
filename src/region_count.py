import time
from numpy import zeros, uint32, uint64, save, load

class RegionCounts(object):
    """records sequence read counts for a genomic region"""
    def __init__(self, length):
        super(RegionCounts, self).__init__()
        
        self.counts = zeros(length, uint32)
    
    def addRead(self, start, end):
        """add counts for range start, end"""
        
        self.counts[start: end] += 1
    
    def save(self, filename, ordered_gene_coords):
        """saves as in numpy binary format"""
        save(filename, self.counts)
    

