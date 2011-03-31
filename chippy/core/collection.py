"""collections handle storage, loading, grouping of data"""

import gzip
import numpy

from chippy.util.util import make_even_groups

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"
__version__ = '0.1'


class RegionCollection(object):
    """store counts, ranks, labels, run arguments from a read count session"""
    def __init__(self, counts=None, ranks=None, labels=None, info=None,
            filename=None):
        super(RegionCollection, self).__init__()
        
        if counts is not None:
            counts = numpy.array(counts)
        
        if ranks is not None:
            ranks = numpy.array(ranks)
            assert ranks.shape[0] == counts.shape[0],\
                'inconsistent data attributes'
        
        if labels is not None:
            labels = numpy.array(labels)
            assert labels.shape[0] == counts.shape[0],\
                'inconsistent data attributes'
        
        self.counts = counts
        self.ranks = ranks
        self.labels = labels
        self.info = info
        
        if filename is not None:
            assert min(counts, ranks, labels, info) is None,\
                            "Conflicting arguments"
            
            self._load(filename)
        
    
    def __str__(self):
        v = 'RegionCollection(num_records=%s; has_ranks=%s; has_labels=%s)'\
                % (len(self.counts), self.ranks is not None,
                                    self.labels is not None)
        return v
    
    def writeToFile(self, filename):
        """writes a gzipped .npy formatted data store"""
        outfile = gzip.GzipFile(filename, 'w')
        save_data = dict(counts=self.counts, ranks=self.ranks,
                         labels=self.labels, info=self.info)
        numpy.save(outfile, save_data)
        outfile.close()
    
    def _load(self, filename):
        """loads attributes from a gzipped, .npy data structure"""
        infile = gzip.GzipFile(filename, 'r')
        data = numpy.load(infile)
        infile.close()
        
        # remember numpy.load() returns and array object
        # numpy.load().tolist() returns a dict ... wtf !
        
        data = data.tolist()
        for name in data:
            self.__dict__[name] = data[name]
        
    
    def _get_keep_indices(self, filtered=None):
        keep = range(self.counts.shape[0])
        if filtered is not None:
            keep = []
            for i in range(self.counts.shape[0]):
                if filtered(self.counts[i]):
                    keep.append(i)
            if len(keep) == 0:
                raise RuntimeError('Filter operation excluded all data!!')
        
        return keep
    
    def getGrouped(self, group_size, filtered=None):
        """Arguments:
            - filtered: a callback function that takes individual counts
              records and returns True if they're to be included, False for
              exclusion
        """
        keep = self._get_keep_indices(filtered=filtered)
        counts = self.counts.take(keep, axis=0)
        
        if self.ranks is not None:
            ranks = self.ranks.take(keep, axis=0)
        else:
            ranks = None
        
        if self.labels is not None:
            labels = self.labels.take(keep, axis=0)
        else:
            labels = None
        
        return counts, ranks, labels
    
    def itergroups(self, group_size, filtered=None):
        """Arguments:
            - filtered: a callback function that takes individual counts
              records and returns True if they're to be included, False for
              exclusion
        """
        counts, ranks, labels = self.getGrouped(group_size, filtered=filtered)
        for i in range(counts.shape[0]):
            if ranks is None:
                rank_data = None
            else:
                rank_data = ranks[i]
            
            if labels is None:
                label_data = None
            else:
                label_data = labels[i]
            
            count_data = counts[i]
            yield count_data, rank_data, label_data
        
    



