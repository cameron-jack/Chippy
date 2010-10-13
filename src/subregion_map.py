from __future__ import division
import re
from glob import glob1
import numpy as np
from os import path as p

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


class MapScores(object):
    """Once the mappability scores have been created, this class allows you to
    load the data to make it useable"""

    def __init__(self, path):
        if p.exists(path):
            mapscores_fns = glob1(path, '*chr?*mapscore*.npy')
            if len(mapscores_fns) != 0:
                self.path = path
                self.mapscores_fns = mapscores_fns
                self.num_files = len(mapscores_fns)
            else:
                raise IOError('Specified path has no valid npy files.')
        else:
            raise IOError('Specified path does not exist.')

        # Create dictionaries (chrom: filename, chrom: map data)
        self.mapscore_dict = {}
        self._chrom_path = {}
        for fn in self.mapscores_fns:
            chrom = re.findall('chr[0-9][0-9]|chr[0-9XY]', fn)[0].strip('chr')
            self._chrom_path[chrom] = p.join(self.path,fn)

    def __str__(self):
        return '%d map-score files found at location: %s' % \
               (self.num_files, self.path)

    def __getitem__(self, chrom):
        chrom = str(chrom)
        try:
            mapscores = self.mapscore_dict[chrom]
        except KeyError:
            self.mapscore_dict[chrom] = np.load(self._chrom_path[chrom])
            mapscores = self.mapscore_dict[chrom]
        return mapscores

    def __delitem__(self, chrom):
        chrom = str(chrom)
        try:
            del(self.mapscore_dict[chrom])
        except KeyError:
            pass

