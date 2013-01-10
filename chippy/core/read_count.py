from __future__ import division
import sys
sys.path.extend(['..'])
from os import path, listdir

import warnings
warnings.filterwarnings('ignore', 'Not using MPI as mpi4py not found')
import numpy, subprocess

from cogent import LoadTable
from cogent.util.progress_display import display_wrap

from chippy.util.definition import NULL_STRAND, PLUS_STRAND, MINUS_STRAND
from chippy.parse.bed import BedRep

__author__ = "Anuj Pahwa, Gavin Huttley, Cameron Jack"
__copyright__ = "Copyright 2012, Gavin Huttley, Cameron Jack, Anuj Pahwa"
__credits__ = ['Anuj Pahwa', 'Gavin Huttley', 'Cameron Jack']
__license__ = "GPL"
__maintainer__ = "Cameron Jack"
__email__ = "cameron.jack@anu.edu.au"
__status__ = "development"
__version__ = '0.1'

@display_wrap
def make_contig_counts(mapped_read_path, max_read_length=None,
                count_max_length=False, strand=NULL_STRAND, sep='\t',
                is_sorted=True, ui=None):
    """returns a numpy array representing read counts
    
    Arguments:
        - mapped_read_path: path to table containing read coordinates,
          frequency data, or the data itself
        - max_read_length: maximum length of a read length
        - count_max_length: if max_read_length provided, all mapped seqs set
          to this length
        - strand: only reads from specified strand are added. Default is both.
        - sep: the delimiter in the read coordinates file
        - is_sorted: whether the read file is already sorted
    """

    if type(mapped_read_path) == str or type(mapped_read_path) == unicode:
        data = LoadTable(mapped_read_path, sep=sep)
        assert list(data.Header) == ['start', 'length', 'strand', 'freq'],\
                "mapped read Table header doesn't match expected"
        if not is_sorted:
            data = data.sorted(columns='start')
        data = data.array.astype(numpy.int32)
    else:
        # It came from a BED file
        data = mapped_read_path

    #print mapped_read_path

    if count_max_length:
        assert max_read_length, 'must specify max_read_length to use'\
                                ' count_max_length'

    total = data.shape[0]
    max_read_length = max_read_length or numpy.inf
    total_length = data[-1][0] + data[-1][1]
    counts = numpy.zeros(total_length, numpy.int32)
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

        if type(mapped_read_path) == str or type(mapped_read_path) == unicode:
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
                    work = numpy.zeros(span, dtype=result.dtype)
                    if slice.step == -1:
                        # if minus strand, we reverse it first
                        work[diff:] = result
                    else:
                        work[:result.shape[0]] = result
                    result = work
                    warnings.warn(msg, UserWarning, 2)
                    
                    # TODO pad the result with zeros?
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
        new = numpy.zeros(length, dtype=self.counts.dtype)
        new[:self.counts.shape[0]] += self.counts
        new[:other.counts.shape[0]] += other.counts
        if self.strand == other.strand:
            strand = self.strand
        else:
            strand = NULL_STRAND
        
        return self.__class__(counts=new, strand=strand)

def read_all_beds(counts_dir):
    """ Read all BED files into memory. We can access the chroms from them late
        via get_chrom_as_nparray(chrom) """

    if type(counts_dir) == str:
        counts_dirs = [counts_dir]
    else:
        counts_dirs = counts_dir

    bed_reps = []
    for counts_dir in counts_dirs:
        dir_list = listdir(counts_dir)
        for file_name in dir_list:
            if file_name[-4:] == '.bed' or file_name[-7:] == '.bed.gz':
                bed_reps.append(BedRep(path.join(counts_dir, file_name)))

    return bed_reps

def get_combined_counts(counts_dir, bed_reps, chrom_name, max_read_length, count_max_length):
    """returns a single WholeChrom, potentially summing counts from different
    sequencer runs"""
    if type(counts_dir) == str:
        counts_dirs = [counts_dir]
    else:
        counts_dirs = counts_dir
    
    chrom = None

    # first check for chrom files and read in
    for counts_dir in counts_dirs:

        chrom_counts_path = path.join(counts_dir,
                    'chr%s.txt.gz' % chrom_name)
        if path.exists(chrom_counts_path):
            print '\t%s' % chrom_counts_path
            counts = WholeChrom(chrom_counts_path,
                            max_read_length=max_read_length,
                            count_max_length=count_max_length)
            if chrom is None:
                chrom = counts
            else:
                chrom = counts + chrom

    # now add data from any BEDs
    for bed_rep in bed_reps:
        chrom_array = bed_rep.get_chrom_as_nparray(chrom_name)
        print '\t%s: chr%s' % (bed_rep.name, chrom_name)
        counts = WholeChrom(chrom_array, max_read_length=max_read_length,
                count_max_length=count_max_length)
        if counts != 0:
            if chrom is None:
                chrom = counts
            else:
                chrom = counts + chrom

    # Now combine counts from all BEDs
    return chrom

def run_command(command):
    """executes a command"""
    PIPE = subprocess.PIPE
    r = subprocess.Popen(command, shell=True, universal_newlines=True,
        stdout=PIPE, stderr=PIPE, bufsize=-1)

    stdout, stderr = r.communicate()
    returncode = r.returncode
    return returncode, stdout, stderr

def readBAM(bamfile_path, ROIs):
    """ get sequence reads for each ROI and add the counts in.
    Relies on Samtools view """
    valid_flags = set([0, 16, 83, 99, 147, 163])
    second_read_flags = set([83,147])
    for roi in ROIs:
        chrom = roi.chrom
        start = roi.start
        end = roi.end
        command = 'samtools view ' + bamfile_path + ' ' + chrom +\
                ':' + start + '-' + end
        returncode, stdout, stderr = run_command(command)
        if returncode == 0:
            bam_lines = stdout.split('\n')
            for record in bam_lines:
                record_parts = record.split('\t')
                record_flags = int(record_parts[1])
                if record_flags in valid_flags:
                    record_start = int(record_parts[3])
                    record_length = len(record_parts[7])
                    if record_flags in second_read_flags:
                        record_strand = MINUS_STRAND
                        record_first = record_start - record_length
                        record_last = record_start
                    else:
                        record_strand = PLUS_STRAND
                        record_first = record_start
                        record_last = record_start + record_length
                    offset_left = record_first - roi.start
                    if offset_left < 0:
                        offset_left = 0
                    offset_right = roi.end - record_last
                    if offset_right > len(roi.count_array):
                        offset_right = len(roi.count_array) - 1
                    roi.count_array[offset_left:offset_right] += 1
    return ROIs


