from __future__ import division

import sys
sys.path.extend(['..'])
import numpy
import uuid

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
            self.counts = numpy.zeros(end - start, dtype=numpy.float32)
        except ValueError:
            print self.chrom, self.label, self.rank, self.TSS,\
                    start, end, self.strand

        # create a unique ID for the ROI
        id_parts = [self.chrom, self.start, self.end, self.strand,
                    str(uuid.uuid4())]
        self.unique_id = '_'.join(map(str, id_parts))

    def __repr__(self):
        return tuple(self.chrom, self.start, self.end)