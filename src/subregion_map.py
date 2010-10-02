from __future__ import division
import numpy as np

class SubregionMap(object):
    """A class that handles the mappability score of a set of coordinates
    based on counts produced by genomic control data"""

    def __init__(self, count_fn, chrom, chrom_length, sites, control_seq_length):
        counts = np.load(count_fn)
        shape = counts.shape
        self.window_size = shape[1]/2 # window size is half the length
        assert (len(sites) == shape[0])
        self.coords = self._compute_gene_coords(sites, self.window_size,
                                                chrom_length)
        self.chrom = chrom
        self.map_scores = self._calcMappabilityScore(counts, control_seq_length)

    def _calcMappabilityScore(self, counts, control_seq_length):
        """Ideally, control sequences of length n will be mapped and counted
        n times. This method returns the ratio of counts for a nucleotide
        position to the ideal number of times it would be mapped."""

        return counts/control_seq_length

    def _compute_gene_coords(self, sites, window_size, chrom_length):
        coords = []
        for tss in sites:
            start = max([tss - window_size, 0])
            end = min([tss + window_size, chrom_length])
            coords += [(start, end)]
        return coords

    def getMappabilityScore(self):
        return self.map_scores

    def saveMappabilityScore(self, filename):
        """Save mapability scores as a numpy binary"""
        np.save(filename,self.map_scores)


#from segment_count import get_gene_coords
#from cogent import LoadTable

#gene_coords_file = '../data/mouse_gene_coords.txt'
#gene_coords = get_gene_coords(gene_coords_file, 'Y')
#start_sites = [tss for (tss, strand) in gene_coords]
#chrom_lengths_file = '../data/mouse_chrom_lengths_release_58.txt'
#chrom_lengths = LoadTable(chrom_lengths_file, sep='\t')
#chrom_lengths = dict(chrom_lengths.getRawData(['chrom', 'length']))

#print chrom_lengths['Y']
#print start_sites
#exit()
#test = SubregionMap('../data/mouse_aln_chrY-window_10000-chrY.npy', 'Y', chrom_lengths['Y'], start_sites, 75)
#mapScore = test.getMapabilityScore()
#test.saveMapabilityScore('../data/mapscores_Y')

