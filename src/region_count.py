from __future__ import division

from os import path as p
from glob import glob1
import re
from numpy import zeros, uint16, uint32, int32, savez, load, array

from cogent import LoadTable
from cogent.util.progress_display import display_wrap

import util
from definition import NULL_STRAND, PLUS_STRAND, MINUS_STRAND

# def _get_start_index(data, last_start_index, region_start, region_end, max_read_length):
#     """returns the start / stop index"""
#     start_index = -1
#     end_index = data.shape[0]
#     
#     begin = max(region_start - max_read_length, 0)
#     limit = region_end + max_read_length
#     
#     for i in range(last_start_index, end_index):
#         s = data[i]
#         if begin < s < limit and start_index < 0:
#             start_index = max(i-1, 0)
#             break
#     
#     return start_index
# 
# def _count_reads(data, counts, start_index, region_start, region_end, region_strand, max_read_length):
#     """data and counts are numpy arrays"""
#     # find the region containing the data, assumes it's sorted
#     for i in range(start_index, data.shape[0]):
#         read_start, read_length, read_strand, read_freq = data[i]
#         if region_strand != 0 and region_strand != read_strand:
#             continue
#         
#         if read_strand == 1:
#             read_end = read_start + read_length
#         else:
#             read_end = read_start + 1
#             read_start = read_end - read_length
#         
#         if region_end <= read_start or read_end < region_start:
#             continue
#         
#         read_start = max(region_start, read_start) - region_start
#         if read_start > region_end:
#             break
#         
#         read_end = min(region_end, read_end) - region_start
#         counts[read_start: read_end] += read_freq
# 
# 
# try:
#     # raise ImportError
#     from _region_count import count_reads, get_start_index
# except ImportError:
#     count_reads = _count_reads
#     get_start_index = _get_start_index
# 
# def CountsForRegion(counts, start, end):
#     """factory function returning a numpy array of counts"""
#     def call((read_start, length, strand, freq)):
#         """returns a numpy array"""
#         if strand == 1:
#             read_end = read_start + length
#         else:
#             read_end = read_start + 1
#             read_start = read_end - length
#         
#         if end <= read_start or read_end < start:
#             return
#         
#         read_start = max(start, read_start) - start
#         read_end = min(end, read_end) - start
#         counts[read_start: read_end] += freq
#         return
#     
#     return call

class WholeChrom(object):
    def __init__(self, mapped_reads, max_read_length, strand=0, sep='\t', is_sorted=True):
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
        self.total_count = zeros(self.data[-1][0] + max_read_length, int32)
    
    def __setitem__(self, slice, value):
        start, end = min([slice.start, slice.stop]), max([slice.start, slice.stop])
        if slice.step is not None:
            if self.strand != slice.step:
                return
        
        self.total_count[start:end] += value
    
    def __getitem__(self, slice):
        start, end = min([slice.start, slice.stop]), max([slice.start, slice.stop])
        return self.total_count[start:end]
    
    @display_wrap
    def update(self, ui=None):
        total = self.data.shape[0]
        for i in range(total):
            if i % 10 == 0:
                ui.display('Adding reads [%d / %d]' % (i, total), i / total)
            
            start, length, strand, freq = self.data[i]
            if self.strand != NULL_STRAND and self.strand != strand:
                continue
            
            if strand == MINUS_STRAND:
                start = start + 1 - length
            end = start + length
            self[start:end] = freq
    
# 
# class MinimalRegionCount(object):
#     """returns read counts for provided coordinates from a delimited file"""
#     def __init__(self, mapped_reads, max_read_length, sep='\t', is_sorted=True):
#         super(MinimalRegionCount, self).__init__()
#         self.max_read_length = max_read_length
#         
#         data = LoadTable(mapped_reads, sep=sep)
#         assert list(data.Header) == ['start', 'length', 'strand', 'freq'],\
#             "mapped read Table header doesn't match expected"
#         
#         if not is_sorted:
#             data = data.sorted(columns='start')
#         
#         self.data = data.array.astype(int32)
#         self.last_start_index = 0
#     
#     def __getitem__(self, slice):
#         """returns counts for region"""
#         start, end = min([slice.start, slice.stop]), max([slice.start, slice.stop])
#         if slice.step is None:
#             strand = 0
#         else:
#             strand = slice.step
#         
#         counts = zeros(end-start, int32)
#         
#         start_index = get_start_index(self.data[:, 0], self.last_start_index,
#                 start, end, self.max_read_length)
#         
#         if start_index < 0:
#             return counts
#         
#         self.last_start_index = start_index
#         
#         count_reads(self.data, counts, start_index, start, end, strand,
#                     self.max_read_length)
#         
#         return counts
#     

class RegionCounts(object):
    """records sequence read counts for a genomic region"""
    def __init__(self, length, one_based=True):
        super(RegionCounts, self).__init__()
        # using uint16 allows for a maximum count value of 65535 which
        # *should* be enough.
        self.length = length
        self.counts = zeros(length, uint16)
        
        # if the mapping software counts from one
        self._adjust = [0, -1][one_based]

    def addRead(self, start, end):
        """add counts for range start, end, adjusting for whether the number
        system starts at 1"""
        self.counts[self._adjust + start: end] += 1

    def _get_subregion_counts(self, ordered_coords, window_size, control=False):
        annotated_counts = zeros((len(ordered_coords), 2 * window_size),
                                 self.counts.dtype)

        for row_index, (tss, strand) in enumerate(ordered_coords):
            if control:
                strand = 1
            stride = strand
            # positive strand
            if strand == 1:
                start = max([tss - window_size, 0])
                end = min([tss + window_size, self.length])
            # negative strand
            else:
                # start which actually be > end
                start = min([tss + window_size, self.length])
                end = max([tss - window_size, 0])

            try:
                annotated_counts[row_index] = self.counts[start: end: stride]
            except ValueError:
                temp = self.counts[start: end: stride].copy()
                temp.resize(2*window_size)
                annotated_counts[row_index] = temp

        return annotated_counts

    def save(self, filename, ordered_coords, window_size, control=False):
        """saves in numpy binary format

        Arguments:
            - filename: full path to be saved to, NOTE: numpy adds the .npy
              suffix
            - ordered_coords: a series with [(tss, strand), ..] where
              tss stands for transcription start site
            - window_size: counts in a window  of tss +/- window_size
              saved. Reverse strand counts are reversed.
            - control: When control is true we are always
              going to assume that the control sequences are on the positive
              strand.
        """
        
        annotated_counts = self._get_subregion_counts(ordered_coords,
                                                      window_size, control)
        
        annotated_counts = zip(map(str, range(len(annotated_counts))), 
                                annotated_counts)
        savez(filename, **dict(annotated_counts))


class IrregularRegionCounts(RegionCounts):
    """returns counts for unequal sized regions"""
    def _get_subregion_counts(self, ordered_coords, control=False):
        annotated_counts = []
        for x, y, strand in ordered_coords:
            stride = strand
            # positive strand
            if strand == 1:
                start = min([x, y])
                end = max([x, y])
            # negative strand
            else:
                start = max([x, y])
                end = min([x, y]) - 1
            
            annotated_counts.append(self.counts[start: end: stride])
        
        return annotated_counts
    
    def save(self, filename, ordered_coords, control=False):
        """saves in numpy binary format
        
        Arguments:
            - filename: full path to be saved to, NOTE: numpy adds the .npy
              suffix
            - ordered_coords: a series with [(start, end, strand), ..]
            - control: When control is true we are always
              going to assume that the control sequences are on the positive
              strand.
        """
        
        annotated_counts = self._get_subregion_counts(ordered_coords, control)
        annotated_counts = zip(map(str, range(len(annotated_counts))), 
                                annotated_counts)
        
        savez(filename, **dict(annotated_counts))
    


class RegularRegionCounts(RegionCounts):
    """returns counts for equivalently sized regions"""
    def _get_subregion_counts(self, ordered_coords, window_size, control=False):
        annotated_counts = zeros((len(ordered_coords), 2 * window_size),
                                 self.counts.dtype)

        for row_index, (tss, strand) in enumerate(ordered_coords):
            if control:
                strand = 1
            stride = strand
            # positive strand
            if strand == 1:
                start = max([tss - window_size, 0])
                end = min([tss + window_size, self.length])
            # negative strand
            else:
                # start which actually be > end
                start = min([tss + window_size, self.length])
                end = max([tss - window_size, 0])

            try:
                annotated_counts[row_index] = self.counts[start: end: stride]
            except ValueError:
                temp = self.counts[start: end: stride].copy()
                temp.resize(2*window_size)
                annotated_counts[row_index] = temp

        return annotated_counts

class _CacheCounts(object):
    """docstring for _CacheCounts"""
    
    def __getitem__(self, chrom):
        chrom = str(chrom)
        try:
            counts = self.count_dict[chrom]
        except KeyError:
            counts = load(self._chrom_path[chrom])
            if hasattr(counts, 'files'):
                # assume just string (int) naming
                names = map(str, range(len(counts.files)))
                data = [counts[name] for name in names]
                try:
                    counts = array(data)
                except ValueError:
                    counts = data # uequal length component arrays
            
            if self.window_size is not None:
                start, end = util.get_centred_coords(counts.shape[1],
                    self.window_size)
                counts = counts[:, start:end]
            
            self.count_dict[chrom] = counts
        return counts
    
    def __delitem__(self, chrom):
        chrom = str(chrom)
        try:
            del(self.count_dict[chrom])
        except KeyError:
            pass
        

class CacheLaneCounts(_CacheCounts):
    """Abstracts the handling of stored RegularRegionCounts"""

    def __init__(self, lane, path, window_size=None):
        super(CacheLaneCounts, self).__init__()

        # The path where the .npy files can be found for this particular lane.
        # Validate that path exists and contains valid npy files.
        if p.exists(path):
            count_filenames = glob1(path, '*s_%s*chr?*.np*'% str(lane))
            if len(count_filenames) != 0:
                self.path = path
                self.count_filenames = count_filenames
                self.num_files = len(count_filenames)
            else:
                raise IOError('Specified path has no valid npz files.')
        else:
            raise IOError('Specified path does not exist.')

        # The data from this lane was used to compute the count statistics.
        self.lane = lane

        # Create a dictionary (chrom: open count file)
        self.count_dict = {}
        self._chrom_path = {}
        for fn in self.count_filenames:
            chrom = re.findall('chr[0-9][0-9]|chr[0-9XY]', fn)[0].strip('chr')
            self._chrom_path[chrom] = p.join(self.path, fn)
        
        self.window_size = window_size
    
    def __str__(self):
        return '%d count files found for lane %s at location: %s' % \
               (self.num_files, str(self.lane), self.path)

    
