import time
from numpy import zeros, uint16, save, load

class RegionCounts(object):
    """records sequence read counts for a genomic region"""
    def __init__(self, length):
        super(RegionCounts, self).__init__()

        # using uint16 allows for a maximum count value of 65535 which *should*
        # be enough.
        self.counts = zeros(length, uint16)

    def addRead(self, start, end):
        """add counts for range start, end"""

        self.counts[start: end] += 1

    def save(self, filename, ordered_gene_coords, windowSize):
        """saves as in numpy binary format"""

        annotatedCounts = zeros((len(ordered_gene_coords), windowSize))

        row_index = 0
        for [chrom, tss, five, three, strand] in ordered_gene_coords:
            # positive strand
            if strand == 1:
                annotatedCounts[row_index] = self.counts[five:three + 1]
            # negative strand
            else:
                annotatedCounts[row_index] = self.counts[three:five + 1]
            row_index += 1

        print 'Saving ' + filename + '\n\n'
        save(filename, annotatedCounts)


