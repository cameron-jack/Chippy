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

def _make_mean(axis):
    def call(data):
        return data.mean(axis=axis)
    return call

def _make_std(axis):
    def call(data):
        return data.std(axis=axis, ddof=1)
    return call

column_mean = _make_mean(0)
column_stdev = _make_std(0)
row_mean = _make_mean(1)
rank_mean = _make_mean(None)

def _get_keep_indices(data, filtered=None):
    keep = range(data.shape[0])
    
    if filtered is not None:
        keep = []
        for i in range(data.shape[0]):
            if filtered(data[i]):
                keep.append(i)
        if len(keep) == 0:
            raise RuntimeError('Filter operation excluded all data!!')
    
    return keep

class _GenericCollection(object):
    def __init__(self, counts=None, ranks=None, labels=None, info=None):
        super(_GenericCollection, self).__init__()
        
        if counts is not None:
            counts = numpy.array(counts)
        
        if ranks is not None:
            ranks = numpy.array(ranks).astype(float)
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
    

class RegionCollection(_GenericCollection):
    """store counts, ranks, labels, run arguments from a read count session"""
    def __init__(self, filename=None, **kwargs):
        super(RegionCollection, self).__init__(**kwargs)
        
        if filename is not None:
            assert min(self.counts, self.ranks, self.labels,
                self.info) is None, "Conflicting arguments"
            
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
            value = data[name]
            self.__dict__[name] = value
            if name == 'ranks' and value is not None:
                self.__dict__[name] = value.astype(float)
        
    
    def normalisedCounts(self, axis=None):
        """Arguments:
            - axis: normalisation is done per column (axis=0), per rows
              (axis=1) or across the entire collection (axis=None).
        """
        copied = self.counts.copy()
        copied = copied.astype(float)
        means = copied.mean(axis=axis)
        stdevs = copied.std(axis=axis, ddof=1)
        if axis == 1:
            means = numpy.vstack(means)
            stdevs = numpy.vstack(stdevs)
        copied -= means
        copied /= stdevs
        return copied
    
    def filteredNormalised(self, cutoff=3.0, axis=None):
        """returns a new RegionCollection excluding records above cutoff"""
        data = self.normalisedCounts(axis=axis)
        func = lambda x: (x < cutoff).all()
        
        indices = _get_keep_indices(data, filtered=func)
        counts = self.counts.take(indices, axis=0)
        if self.ranks is not None:
            ranks = self.ranks.take(indices, axis=0)
        else:
            ranks = None
        
        if self.labels is not None:
            labels = self.labels.take(indices, axis=0)
        else:
            labels = None
        
        if self.info is None:
            info = None
        else:
            info = self.info.copy()
            info['filteredNormalised'] = cutoff
        
        new = self.__class__(counts=counts, labels=labels, ranks=ranks,
                    info=info)
        return new
    
    def getGrouped(self, group_size, filtered=None, normalised=False,
                                axis=None, indices=None):
        """Arguments:
            - filtered: a callback function that takes individual counts
              records and returns True if they're to be included, False for
              exclusion
            - normalised: if True, the results of normalisedCounts are
              supplied to filtered
            - axis: normalisation is done per column (axis=0), per rows
              (axis=1) or across the entire collection (axis=None).
            - indices: indices of rows to keep
        """
        
        assert not (normalised and indices), \
                    'Normalisation only relevant if indices not provided'
        
        if normalised:
            data = self.normalisedCounts(axis=axis)
        else:
            data = self.counts
        
        if not indices:
            indices = _get_keep_indices(data, filtered=filtered)
        
        counts = self.counts.take(indices, axis=0)
        
        if group_size > 1:
            counts = make_even_groups(counts, group_size)
        
        if self.ranks is not None:
            ranks = self.ranks.take(indices, axis=0)
            ranks = make_even_groups(ranks, group_size)
        else:
            ranks = None
        
        if self.labels is not None:
            labels = self.labels.take(indices, axis=0)
            labels = make_even_groups(labels, group_size)
        else:
            labels = None
        
        return numpy.array(counts), numpy.array(ranks), numpy.array(labels)
    
    def itergroups(self, **kwargs):
        """Arguments same as for getGrouped
        """
        counts, ranks, labels = self.getGrouped(**kwargs)
        for i in range(len(counts)):
            if ranks is None:
                rank_data = None
            else:
                rank_data = numpy.array(ranks[i])
            
            if labels is None:
                label_data = None
            else:
                label_data = numpy.array(labels[i])
            
            count_data = numpy.array(counts[i])
            yield count_data, rank_data, label_data
        
    
    def iterDescriptiveStats(self, **kwargs):
        """return column means & stdevs for each group"""
        for counts, ranks, labels in self.itergroups(**kwargs):
            yield column_mean(counts), column_stdev(counts), rank_mean(ranks)
    
    def iterTransformedGroups(self, rank_func=rank_mean,
                    counts_func=column_mean, **kwargs):
        """docstring for transformedGroups"""
        for counts, ranks, labels in self.itergroups(**kwargs):
            c = counts_func(counts)
            r = rank_func(ranks)
            yield c, r
        
    def filteredByLabel(self, labels):
        """returns a new collection object consisting of just the subset"""
        if self.labels is None:
            raise RuntimeError('No labels')
        
        if type(labels) == str:
            labels = [labels]
        
        if self.ranks is None:
            rank_data = None
        else:
            rank_data = []
        
        label_data = []
        count_data = []
        counts = self.counts
        for i in range(counts.shape[0]):
            if self.labels[i] in labels:
                count_data.append(self.counts[i])
                label_data.append(self.labels[i])
                if rank_data is not None:
                    rank_data.append(self.ranks[i])
        
        return self.__class__(counts=count_data, ranks=rank_data,
                            labels=label_data, info=self.info)
    

