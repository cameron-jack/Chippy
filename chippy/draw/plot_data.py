from __future__ import division

import sys, numpy
sys.path.extend(['..', '../src'])
from math import log
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

class PlotPoint(object):
    """ A simple object representation of a point for plotting """
    def __init__(self, x=0, y=0):
        super(PlotPoint, self).__init__()
        self.x = x
        self.y = y

    def __repr__(self):
        return repr((self.x, self.y))

    def get_logX(self):
        """ safe log base 2 transform self.x """
        if self.x <= 0:
            x = 0
        else:
            x = log(self.x + 1, 2)
        return x

    def get_logY(self):
        """ safe log base 2 transform self.y """
        if self.y <= 0:
            y = 0
        else:
            y = log(self.y + 1, 2)
        return y

class PlotLine(object):
    """ Represents an individual line to be plotted.
        Can be smoothed, binned and sorted.
        Comprised of data to be plotted, rank order of data, a label
        and the name of the originating study.
        Optional data members are color and stderr
    """
    def __init__(self, counts, rank=None, label=None, study=None, **kwargs):
        super(PlotLine, self).__init__()
        self.counts = counts
        self.rank = rank
        if type(label) == list or type(label) == tuple:
            label = ', '.join(label)
        self.label = label
        self.study = study
        if 'color' in kwargs:
            self.color = kwargs['color']
        else:
            self.color = 'grey'
        if 'alpha' in kwargs:
            self.alpha = kwargs['alpha']
        else:
            self.alpha = 0.9
        if 'stderr' in kwargs:
            self.stderr = kwargs['stderr']
        else:
            self.stderr = None

    def __repr__(self):
        return ', '.join(map(str,
                [self.counts, self.rank, self.label, self.study]))

    def applySmoothing(self, smooth_width):
        """ User-directed binning or smoothing of count line
        data is done here """
        self.counts = smooth(self.counts, smooth_width)

    def applyBinning(self, bin_width):
        """ For every bin_width, sum the counts. Output array size is
            same, just filled with mean values of each bin - giving a
            normalised score per base. Bin_width must be an integer
            factor of the window size.
        """
        rr = RunRecord('apply_binning')
        if bin_width and bin_width > 0:
            if len(self.counts)%bin_width:
                rr.dieOnCritical('Bin width is not an integer '+\
                        'factor of window size', bin_width)
            tmp_array = numpy.array(self.counts)
            for k in range(0, len(self.counts), bin_width):
                bin_sum = 0
                for i in xrange(bin_width):
                    bin_sum += self.counts[k+i]
                for i in xrange(bin_width):
                    tmp_array[k+i] = bin_sum/bin_width
            self.counts = tmp_array

    def countsAsArray(self):
        """ Return self as numpy array of Y coords """
        if type(self.counts) == 'numpy.ndarray':
            return self.counts
        else:
            return numpy.array(self.counts)

    def getMaxCount(self, include_stderr=False, se_adjust=1):
        if len(self.counts) > 0:
            if not include_stderr or self.stderr is None:
                return max(self.counts)
            else:
                return max(self.counts) + self.stderr * se_adjust
        return None

    def getMinCount(self, include_stderr=False, se_adjust=1):
        if len(self.counts) > 0:
            if not include_stderr or self.stderr is None:
                return min(self.counts)
            else:
                return min(self.counts) - self.stderr * se_adjust
        return None

    def getLabelsAsList(self):
        if type(self.label) == str:
            labels = self.label.split(', ')
        elif type(self.label) == numpy.ndarray:
            labels = self.label.tolist()

        return labels

