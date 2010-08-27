from numpy import zeros, uint16, save, load

class RegionCounts(object):
    """records sequence read counts for a genomic region"""
    def __init__(self, length, one_based=True):
        super(RegionCounts, self).__init__()
        # using uint16 allows for a maximum count value of 65535 which *should*
        # be enough.
        self.length = length
        self.counts = zeros(length, uint16)
        
        # if the mapping software counts from one
        self._adjust = [0, -1][one_based]

    def addRead(self, start, end):
        """add counts for range start, end, adjusting for whether the number
        system starts at 1"""
        self.counts[self._adjust + start: end] += 1
    
    def _get_subregion_counts(self, ordered_coords, window_size):
        annotated_counts = zeros((len(ordered_coords), 2 * window_size), uint16)
        
        for row_index, (tss, strand) in enumerate(ordered_coords):
            stride = strand
            # positive strand
            if strand == 1:
                start = max([tss - window_size, 0])
                end = min([tss + window_size, self.length])
            # negative strand
            else:
                # start which actually be > end
                start = min([tss + window_size, self.length])
                end = max([tss - window_size, 0])
            annotated_counts[row_index] = self.counts[start: end: stride]
        return annotated_counts
    
    def save(self, filename, ordered_coords, window_size):
        """saves in numpy binary format
        
        Arguments:
            - filename: full path to be saved to, NOTE: numpy adds the .npy
              suffix
            - ordered_coords: a series with [(tss, strand), ..] where
              tss stands for transcription start site
            - window_size: counts in a window  of tss +/- window_size
              saved. Reverse strand counts are reversed.
        """
        
        annotated_counts = self._get_subregion_counts(ordered_coords,
                                                      window_size)
        save(filename, annotated_counts)

