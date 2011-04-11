from __future__ import division

from os import path as p
from glob import glob1
import re
import warnings
warnings.filterwarnings('ignore', 'Not using MPI as mpi4py not found')

from numpy import zeros, uint16, uint32, int32, savez, load, array, inf

from cogent import LoadTable
from cogent.util.progress_display import display_wrap

from chippy.util import util
from chippy.util.definition import NULL_STRAND, PLUS_STRAND, MINUS_STRAND

__author__ = "Anuj Pahwa, Gavin Huttley"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley"
__credits__ = ["Anuj Pahwa", "Gavin Huttley"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"
__version__ = '0.1'

class WholeChrom(object):
    def __init__(self, mapped_reads, max_read_length=None, strand=0, sep='\t', is_sorted=True):
        super(WholeChrom, self).__init__()
        self.max_read_length = max_read_length
        
        data = LoadTable(mapped_reads, sep=sep)
        assert list(data.Header) == ['start', 'length', 'strand', 'freq'],\
            "mapped read Table header doesn't match expected"
        
        if not is_sorted:
            data = data.sorted(columns='start')
        
        self.data = data.array.astype(int32)
        self.last_start_index = 0
        self.strand = strand
        
        total_length = self.data[-1][0] + self.data[-1][1]
        self.total_count = zeros(total_length, int32)
    
    def __setitem__(self, slice, value):
        # this is to handle the construction of a strand specific count contig
        start, end = min([slice.start, slice.stop]), max([slice.start, slice.stop])
        if slice.step is not None:
            if self.strand != slice.step:
                # if we defined a strand for this instance, we only add reads
                # that have the same strand
                return
        try:
            value.shape[0]
            return # numpy modifies arrays in place
        except IndexError:
            pass
        
        self.total_count[start:end] += value
    
    def __getitem__(self, slice):
        # not clear how to handle a slice that returns an shorter array
        msg = "you've sliced beyond the limits of the contig"
        try:
            result = self.total_count[slice]
            if slice.stop is not None:
                if abs(slice.stop - slice.start) != result.shape[0]:
                    warnings.warn(msg, UserWarning, 2)
        except IndexError:
            warnings.warn(msg, UserWarning, 2)
            result = None
        
        return result
    
    @display_wrap
    def update(self, ui=None):
        """applies referenced read count data to produce a single numpy array"""
        total = self.data.shape[0]
        max_read_length = self.max_read_length or inf
        for i in range(total):
            if i % 10 == 0:
                ui.display('Adding reads [%d / %d]' % (i, total), i / total)
            
            start, length, strand, freq = self.data[i]
            
            if self.strand != NULL_STRAND and self.strand != strand:
                continue
            
            end = start + length
            if max_read_length < length:
                diff = length - max_read_length
                if strand == PLUS_STRAND:
                    end -= diff
                elif strand == MINUS_STRAND:
                    start += diff
            
            self[start:end] += freq
        
    
