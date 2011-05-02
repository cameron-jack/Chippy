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

def column_sum(data):
    """returns the column sums"""
    assert len(data.shape) == 2
    return data.sum(axis=0)

def _make_mean(axis):
    def call(data):
        return data.mean(axis=axis)
    return call

def _make_std(axis):
    def call(data):
        return data.std(axis=axis, ddof=1)
    return call

def chebyshev_upper(p):
    """returns k such that the probability that a variable will be greater
    than k standard deviations from its mean is <= p. This is Chebyshev's
    one-sided inequality."""
    return numpy.sqrt(1/p - 1)

def normalised_data(data, axis=None):
    """returns a new normalised array
    
    Arguments:
        - axis: axis on which to compute mean and std dev. If None, all data
          used. If 0/1, row-wise / column-wise statistics are computed.
    """
    copied = data.copy()
    copied = copied.astype(float)
    means = copied.mean(axis=axis)
    stdevs = copied.std(axis=axis, ddof=1)
    if axis == 1:
        means = numpy.vstack(means)
        stdevs = numpy.vstack(stdevs)
    copied -= means
    copied /= stdevs
    return copied

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
        self.N = None
        if counts is not None:
            self.N = self.counts.shape[0]
    
    @property
    def Total(self):
        """returns the total sum of counts"""
        return self.counts.sum()
    

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
        
        self.N = self.counts.shape[0]
    
    def asfloats(self):
        """returns new RegionCollection with counts as floats"""
        counts = self.counts.astype(float)
        return self.__class__(counts=counts, ranks=self.ranks,
                            labels=self.labels, info=self.info)
        
    def asfreqs(self):
        """returns new instance with each count a frequency of Total"""
        new = self.asfloats()
        # now just divide by total
        new.counts /= self.Total
        return new
    
    def normalised(self, axis=None):
        """returns new RegionCollection with counts normalised
        
        Arguments:
            - axis: normalisation is done per column (axis=0), per rows
              (axis=1) or across the entire collection (axis=None).
        """
        counts = normalised_data(self.counts, axis=axis)
        new = self.__class__(counts=counts, labels=self.labels,
                ranks=self.ranks, info=self.info)
        
        return new
    
    def transformed(self, rank_func=rank_mean, counts_func=column_mean):
        """transforms all counts and ranks"""
        c = counts_func(self.counts)
        if self.ranks is not None:
            r = rank_func(self.ranks)
        else:
            r = None
        return c, r
    
    def take(self, indices):
        """returns new instance corresponding to just the indices"""
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
        
        new = self.__class__(counts=counts, labels=labels, ranks=ranks,
                    info=info)
        
        return new
    
    def filtered(self, func, axis=None):
        """returns new RegionCollection excluding records where func is False"""
        indices = _get_keep_indices(self.counts, func)
        return self.take(indices)
    
    def filteredChebyshevUpper(self, p=0.05, axis=None):
        """returns a new RegionCollection excluding records with excessive
        reads using a one-sided Chebyshev's inequality"""
        if not (0 <= p <= 1):
            raise RuntimeError('Probability argument not a valid probability')
        
        k = chebyshev_upper(p)
        if axis is None:
            # only bother computing normalised score for max of each
            # row
            data = self.counts.max(axis=1)
            mean = self.counts.mean()
            stdev = self.counts.std(ddof=1)
            data -= mean
            data /= stdev
            indices = data < k
            data = self.counts[indices]
            if self.labels is not None:
                labels = self.labels[indices]
            else:
                labels = None
            
            if self.ranks is not None:
                ranks = self.ranks[indices]
            else:
                ranks = None
            new = self.__class__(counts=data, ranks=ranks, labels=labels,
                info=self.info)
        else:
            data = normalised_data(self.counts, axis=axis)
            func = lambda x: (x < k).all()
            indices = _get_keep_indices(data, filtered=func)
            new = self.take(indices)
            
        
        if self.info is None:
            info = {'filteredChebyshevUpper': p}
        else:
            info = self.info.copy()
            info['filteredChebyshevUpper'] = p
        
        if new.info:
            new.info.update(info)
        else:
            new.info = info
        
        return new
    
    def getGrouped(self, group_size):
        """returns counts, ranks, labels in group_size"""
        counts = self.counts
        ranks = self.ranks
        labels = self.labels
        if group_size > 1:
            # put into groups
            num_groups, rmdr = divmod(self.N, group_size)
            total = num_groups * group_size
            num_positions = counts.shape[1]
            
            counts = counts[:total]
            counts = counts.reshape(num_groups, group_size, num_positions)
            if ranks is not None:
                ranks = ranks[:total]
                ranks = ranks.reshape(num_groups, group_size)
            
            if labels is not None:
                labels = labels[:total]
                labels = labels.reshape(num_groups, group_size)
        
        counts = numpy.array(counts)
        if ranks is not None:
            ranks = numpy.array(ranks)
        
        if labels is not None:
            labels = numpy.array(labels)
        
        return counts, ranks, labels
    
    def itergroups(self, group_size):
        counts, ranks, labels = self.getGrouped(group_size)
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
        
    
    def iterDescriptiveStats(self, group_size):
        """return column means & stdevs for each group"""
        for counts, ranks, labels in self.itergroups(group_size):
            yield column_mean(counts), column_stdev(counts), rank_mean(ranks)
    
    def iterTransformedGroups(self, group_size, rank_func=rank_mean,
                    counts_func=column_mean):
        for counts, ranks, labels in self.itergroups(group_size):
            c = counts_func(counts)
            r = rank_func(ranks)
            yield c, r, labels
    
    def filteredByLabel(self, labels):
        """returns a new collection object with data corresponding to the
        provided labels"""
        if self.labels is None:
            raise RuntimeError('No labels')
        
        if type(labels) == str:
            labels = [labels]
        
        # determine label indices and use self.take
        indices = []
        for i in range(self.counts.shape[0]):
            if self.labels[i] in labels:
                indices.append(i)
        
        return self.take(indices)
    

