from __future__ import division

import sys, numpy
sys.path.extend(['..', '../src'])

from util import smooth
from chippy.util.run_record import RunRecord

__author__ = 'Cameron Jack'
__copyright__ = 'Copyright 2011-2013, Gavin Huttley, Cameron Jack, Anuj Pahwa'
__credits__ = ['Cameron Jack']
__license__ = 'GPL'
__maintainer__ = 'Cameron Jack'
__email__ = 'cameron.jack@anu.edu.au'
__status__ = 'pre-release'
__version__ = '0.1'

class PlotData(object):
    """ Abstract Base Class: Represents a data entity to be plotted. """
    def __init__(self, *args, **kwargs):
        super(PlotData, self).__init__(*args, **kwargs)
        self.counts = []
        self.rank = 0 # maybe should use -1
        self.label = ''
        self.color = 'grey'

    def countsAsArray(self):
        """ Return self as numpy array of Y coords """
        if type(self.counts) == 'numpy.ndarray':
            return self.counts
        else:
            return numpy.array(self.counts)

    def getMaxCount(self):
        if len(self.counts) > 0:
            return max(self.counts)
        return None

class PlotLine(PlotData):
    """ Represents an individual line to be plotted.
        Can be smoothed, binned and sorted.
        Comprised of data to be plotted, rank order of data, a label
        and the name of the originating study.
        Optional data members are color and stderr
    """
    def __init__(self, counts, rank=None, label=None, study=None,
            *args, **kwargs):
        super(PlotLine, self).__init__(*args, **kwargs)
        self.counts = counts
        self.rank = rank
        self.label = label
        self.study = study
        if 'color' in kwargs:
            self.color = kwargs['color']
        else:
            self.color = 'grey'
        if 'stderr' in kwargs:
            self.stderr = kwargs['stderr']
        else:
            self.stderr = None

    def __repr__(self):
        return self.counts, self.rank, self.label

    def applySmoothing(self, smooth_width):
        """ User-directed binning or smoothing of count line
        data is done here """
        self.counts = smooth(self.counts, smooth_width)

    def applyBinning(self, bin_width):
        """ For every bin_width, sum the counts. Output array size is
            same, just filled with sum values. Bin_width must be an integer
            factor of the window size
        """
        rr = RunRecord('apply_binning')
        if bin_width and bin_width > 0:
            if len(self.counts)%bin_width:
                rr.dieOnCritical('Bin width is not an integer factor of window size', bin_width)
            tmp_array = numpy.array(self.counts)
            for k in self.counts[:-bin_width:bin_width]:
                bin_sum = 0
                for i in xrange(bin_width):
                    bin_sum += self.counts[k+i]
                for i in xrange(bin_width):
                    tmp_array[k+i] = bin_sum
        self.counts = tmp_array

