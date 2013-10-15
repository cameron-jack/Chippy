from __future__ import division

import sys
sys.path.extend(['..'])
import numpy

from chippy.util.definition import NULL_STRAND, PLUS_STRAND, MINUS_STRAND

__author__ = 'Cameron Jack'
__copyright__ = 'Copyright 2011-2013, Anuj Pahwa, Gavin Huttley, Cameron Jack'
__credits__ = ['Cameron Jack']
__license__ = 'GPL'
__maintainer__ = 'Cameron Jack'
__email__ = 'cameron.jack@anu.edu.au'
__status__ = 'Pre-release'
__version__ = '0.1'

class ROI(object):
    """ Maintains and controls data for any genome region for which we
        are collecting data.
    """
    def __init__(self, gene, start, end, label=None):
        """ All coordinates are in Python 0-based [,) space """
        super(ROI, self).__init__()
        #self.species = species
        self.chrom = gene.chrom
        self.gene_id = gene.ensembl_id
        if not label:
            self.label = gene.ensembl_id
        else:
            self.label = label
        self.rank = gene.Rank
        self.TSS = gene.Tss
        self.strand = gene.strand
        self.start = start
        self.end = end
        try:
            self.counts = numpy.zeros(end - start, dtype=numpy.uint32)
        except ValueError:
            print self.chrom, self.label, self.rank, self.TSS,\
                    start, end, self.strand

    def add_counts_to_ROI(self, entry_start, entry_end):
        """ All coordinates are in Python 0-offset space """

        if self.strand == PLUS_STRAND:
            offset_left = entry_start - self.start
            offset_right = entry_end - self.start
        else:
            offset_left = self.end - entry_end
            offset_right = self.end - entry_start

        assert offset_left < offset_right, 'left slicing coord must be '+\
                'less than right slicing coord'

        # Limit the offsets to (0, len(roi.counts))
        offset_left = max(offset_left, 0)
        offset_right = min(offset_right, len(self.counts))

        self.counts[offset_left:offset_right] += 1
        bases_counted = offset_right - offset_left
        return self.counts, bases_counted
