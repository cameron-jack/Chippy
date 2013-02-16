from __future__ import division

from cogent.parse.table import ConvertFields, SeparatorFormatParser
import sys, warnings, re, gzip
import numpy as np

warnings.filterwarnings('ignore', 'Not using MPI as mpi4py not found')
sys.path.extend(['..', '../src'])

__author__ = "Cameron Jack, Gavin Huttley"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley, Cameron Jack"
__credits__ = ["Cameron Jack, Gavin Huttley"]
__license__ = "GPL"
__maintainer__ = "Cameron Jack"
__email__ = "cameron.jack@anu.edu.au"

__status__ = "pre-release"
__version__ = '0.1'

# The following describes a BED file. ChipPy uses BED6.

# The first three required BED fields are: (BED3)
# chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
# chromStart - The starting position of the feature in the chromosome or scaffold.
        # The first base in a chromosome is numbered 0.
# chromEnd - The ending position of the feature in the chromosome or scaffold.
        # The chromEnd base is not included in the display of the feature. For example, the first 100 bases
        # of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.

# The 9 additional optional BED fields are:(BED6)
# name - Defines the name of the BED line. This label is displayed to the left of the BED line in the
        # Genome Browser window when the track is open to full display mode or directly to the left
        # of the item in pack mode.
# score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this
        # annotation data set, the score value will determine the level of gray in which this feature
        # is displayed (higher numbers = darker gray). This table shows the Genome Browser's
        # translation of BED score values into shades of gray:
# strand - Defines the strand - either '+' or '-'.

# (BED12) - these are all group display options and are not used by ChipPy
# thickStart - The starting position at which the feature is drawn thickly (for example,
        # the start codon in gene displays).
# thickEnd - The ending position at which the feature is drawn thickly (for example,
        # the stop codon in gene displays).
# itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb
        # attribute is set to "On", this RBG value will determine the display color of the data
        # contained in this BED line. NOTE: It is recommended that a simple color scheme (eight
        # colors or less) be used with this attribute to avoid overwhelming the color resources
        # of the Genome Browser and your Internet browser.
# blockCount - The number of blocks (exons) in the BED line.
# blockSizes - A comma-separated list of the block sizes. The number of items in this list
        # should correspond to blockCount.
# blockStarts - A comma-separated list of block starts. All of the blockStart positions should
        # be calculated relative to chromStart. The number of items in this list should
        # correspond to blockCount.

def _get_strand(val):
    """ returns 1/-1 for strand in place of '+' or '-' """
    strand = [-1,1][val == '+']
    return strand

pattern = re.compile(r'[0-9,X,Y,MT]+')
def _get_chrom(val):
    """ returns the int component of a chromosome number """
    chrom = pattern.search(val).group(0)
    return chrom

# BED3 defines: chrom, chromStart, chromEnd
bed3_converter = ConvertFields([(0, _get_chrom), (1, int), (2, int)])

# BED6 adds: Name, score, strand
converter = ConvertFields([(0, _get_chrom), (1, int),(2, int), (3, str), (4, int),
                            (5, _get_strand)])

# BED12 additional fields: thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts
complete_converter = ConvertFields([(0, _get_chrom), (1, int),(2, int), (3, str),
        (4, int), (5, _get_strand), (6, int), (7, int), (8, tuple), (9, int),
        (10, tuple), (11, tuple)])

def MinimalBedParser(data, converter=converter):
    """returns data lines from a BED file

    NOTE: BED uses 0-based numbering"""
    # If given a filename for the data
    if type(data) == str:
        if data.endswith('.bed.gz'):
            data = gzip.GzipFile(data, 'rb')
        else:
            data = open(data, 'r')

    header_lines = 0
    data_lines = []

    for row in data:
        if not row.startswith('chr'):
            header_lines += 1
        else:
            data_lines.append(row)


    parser = SeparatorFormatParser(converter=converter, with_header=False,
                    sep="\t")

    for row in parser(data_lines):
        yield row

class BedRep:
    """ defines a BED file. Use a dict of chroms containing a list of tuples.
        This class allows access the same way as for Pickled arrays. The user
        requests a numpy array by chromosome."""

    def __init__(self, bed_file_path):
        parser = MinimalBedParser(bed_file_path)
        self.chrom_dict = {} # only append if not already in list
        self.name = bed_file_path

        for record in parser:
            # Add list of track if chromosome not yet seen
            if record[0] not in self.chrom_dict:
                self.chrom_dict[record[0]] = []

            # add tuple of start, length, strand and frequency of 1, to chrom-specific list
            track = (record[1], record[2]-record[1], record[5], 1)
            self.chrom_dict[record[0]].append(track)

    def get_chrom_as_nparray(self, chrom):
        if chrom not in self.chrom_dict:
            return np.array((0))
        else:
            list_of_chrom_tuples = self.chrom_dict[chrom]
            out_array = np.array(list_of_chrom_tuples)
            return out_array