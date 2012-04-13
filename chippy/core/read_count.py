from __future__ import division

from os import path
from glob import glob1
import re
import warnings
warnings.filterwarnings('ignore', 'Not using MPI as mpi4py not found')

from numpy import zeros, uint16, uint32, int32, savez, load, array, inf

from cogent import LoadTable
from cogent.util.progress_display import display_wrap

from chippy.util import util
from chippy.ref.util import chroms
from chippy.util.definition import NULL_STRAND, PLUS_STRAND, MINUS_STRAND
from chippy.parse.bed import MinimalBedParser

__author__ = "Anuj Pahwa, Gavin Huttley"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley"
__credits__ = ["Anuj Pahwa", "Gavin Huttley"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"
__version__ = '0.1'

@display_wrap
def make_contig_counts(mapped_read_path, max_read_length=None,
                count_max_length=False, strand=NULL_STRAND, sep='\t',
                is_sorted=True, ui=None):
    """returns a numpy array representing read counts
    
    Arguments:
        - mapped_read_path: path to table containing read coordinates,
          frequency data
        - max_read_length: maximum length of a read length
        - count_max_length: if max_read_length provided, all mapped seqs set
          to this length
        - strand: only reads from specified strand are added. Default is both.
        - sep: the delimiter in the read coordinates file
        - is_sorted: whether the read file is already sorted
    """
    data = LoadTable(mapped_read_path, sep=sep)
    assert list(data.Header) == ['start', 'length', 'strand', 'freq'],\
        "mapped read Table header doesn't match expected"
    
    if not is_sorted:
        data = data.sorted(columns='start')
    
    if count_max_length:
        assert max_read_length, 'must specify max_read_length to use'\
                                ' count_max_length'
    data = data.array.astype(int32)
    total = data.shape[0]
    max_read_length = max_read_length or inf
    total_length = data[-1][0] + data[-1][1]
    counts = zeros(total_length, int32)
    for i in range(total):
        if i % 10 == 0:
            ui.display('Adding reads [%d / %d]' % (i, total), i / total)
        
        start, length, read_strand, freq = data[i]
        
        if strand != NULL_STRAND and strand != read_strand:
            continue
        
        end = start + length
        if max_read_length < length:
            diff = length - max_read_length
            if read_strand == PLUS_STRAND:
                end -= diff
            elif read_strand == MINUS_STRAND:
                start += diff
        
        if count_max_length:
            counts[start:start+max_read_length] += freq
        else:
            counts[start:end] += freq
    
    return counts

class WholeChrom(object):
    def __init__(self, mapped_read_path=None, counts=None,
          max_read_length=None, count_max_length=False, strand=NULL_STRAND,
          sep='\t', is_sorted=True):
        super(WholeChrom, self).__init__()
        self.last_start_index = 0
        self.strand = strand
        
        assert not (mapped_read_path and counts),\
                "Cannot provide a counts array AND a path to create one from"
        
        if counts is not None:
            self.counts = counts
        else:
            self.counts = make_contig_counts(mapped_read_path,
                max_read_length=max_read_length,
                count_max_length=count_max_length, strand=strand, sep=sep,
                is_sorted=is_sorted)
        
    
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
            if not (start < self.counts.shape[0] <= end):
                # numpy modifies arrays in place
                return
            else:
                # but in this case we're beyond the end
                diff = self.counts.shape[0] - start
                value = value[:diff]
        except IndexError:
            pass
        
        self.counts[start:end] += value
    
    def __getitem__(self, slice):
        # not clear how to handle a slice that returns an shorter array
        msg = "you've sliced beyond the limits of the contig"
        try:
            result = self.counts[slice]
            if slice.stop is not None:
                span = abs(slice.stop - slice.start)
                if span != result.shape[0]:
                    diff = span - result.shape[0]
                    work = zeros(span, dtype=result.dtype)
                    if slice.step == -1:
                        # if minus strand, we reverse it first
                        work[diff:] = result
                    else:
                        work[:result.shape[0]] = result
                    result = work
                    warnings.warn(msg, UserWarning, 2)
                    
                    # TODO pad the result with zeros
        except IndexError:
            warnings.warn(msg, UserWarning, 2)
            result = None
        
        return result
    
    def __add__(self, other):
        """adds counts returns instance of self.__class__"""
        if not isinstance(other, self.__class__):
            raise RuntimeError('Cannot add %s to %s' % (self.__class__,
                    other.__class__))
        
        length = max(self.counts.shape[0], other.counts.shape[0])
        new = zeros(length, dtype=self.counts.dtype)
        new[:self.counts.shape[0]] += self.counts
        new[:other.counts.shape[0]] += other.counts
        if self.strand == other.strand:
            strand = self.strand
        else:
            strand = NULL_STRAND
        
        return self.__class__(counts=new, strand=strand)
    

def get_combined_counts(counts_dir, chrom_name, max_read_length, count_max_length):
    """returns a single WholeChrom, potentially summing counts from different
    sequencer runs"""
    if type(counts_dir) == str:
        counts_dirs = [counts_dir]
    else:
        counts_dirs = counts_dir
    
    chrom = None
    for counts_dir in counts_dirs:
        chrom_counts_path = path.join(counts_dir,
                    'chr%s.txt.gz' % chrom_name)
        print '\t%s' % chrom_counts_path
        counts = WholeChrom(chrom_counts_path,
                            max_read_length=max_read_length,
                            count_max_length=count_max_length)
        if chrom is None:
            chrom = counts
        else:
            chrom = counts + chrom
    
    return chrom
