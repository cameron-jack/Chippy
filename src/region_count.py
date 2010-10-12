from os import path as p
from glob import glob1
import re
from numpy import zeros, uint16, save, load

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
            - control:When control is true we are always
              going to assume that the control sequences are on the positive
              strand.
        """

        annotated_counts = self._get_subregion_counts(ordered_coords,
                                                      window_size, control)
        save(filename, annotated_counts)


class CacheLaneCounts(object):
    """Abstracts the handling of stored RegionCounts"""

    def __init__(self, lane, path):
        super(CacheLaneCounts, self).__init__()

        # The path where the .npy files can be found for this particular lane.
        # Validate that path exists and contains valid npy files.
        if p.exists(path):
            count_filenames = glob1(path, '*s_%s*chr?*.npy'%str(lane))
            if len(count_filenames) != 0:
                self.path = path
                self.count_filenames = count_filenames
                self.num_files = len(count_filenames)
            else:
                raise IOError('Specified path has no valid npy files.')
        else:
            raise IOError('Specified path does not exist.')

        # The data from this lane was used to compute the count statistics.
        self.lane = lane

        # Create a dictionary (chrom: open count file)
        self.count_dict = {}
        self._chrom_path = {}
        for fn in self.count_filenames:
            chrom = re.findall('chr[0-9][0-9]|chr[0-9XY]', fn)[0].strip('chr')
            self._chrom_path[chrom] = p.join(self.path,fn)

    def __str__(self):
        return '%d count files found for lane %s at location: %s' % \
               (self.num_files, str(self.lane), self.path)

    def __getitem__(self, chrom):
        chrom = str(chrom)
        try:
            counts = self.count_dict[chrom]
        except KeyError:
            self.count_dict[chrom] = load(self._chrom_path[chrom])
            counts = self.count_dict[chrom]
        return counts
    
    def __delitem__(self, chrom):
        chrom = str(chrom)
        try:
            del(self.count_dict[chrom])
        except KeyError:
            pass







