from __future__ import division
import numpy

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"
__version__ = '0.1'


def smoothed(x, num_bins, window_len):
    """returns x and y, with y smoothed"""
    min_x, max_x = min(x), max(x)
    diff = max_x - min_x
    step = diff/num_bins
    h, b = numpy.histogram(x, bins=numpy.arange(min_x, max_x, step), normed=True, new=True)
    x = [(b[i-1]+b[i])/2 for i in range(1, len(b),1)]
    return x, smooth(h, window_len)

def smooth(x, window_len=8, window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal
        window_len: the dimension of the smoothing window
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.
    
    output:
        the smoothed signal
    
    example:
    
    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also:
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
    
    TODO: the window parameter could be the window itself if an array instead of a string
    """
    
    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."
    
    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."
    
    if window_len<3:
        return x
    
    if window not in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', "\
                          "'bartlett', 'blackman'"
    
    s=numpy.r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')
    
    y = numpy.convolve(w/w.sum(),s,mode='same')
    return y[window_len-1:-window_len+1]

