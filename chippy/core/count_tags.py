"""returns tag counts for a specified collection of genes"""
from __future__ import division

import sys
sys.path.extend(['..'])
import numpy

from cogent.util.progress_display import display_wrap
from chippy.core.read_count import get_region_counts
from chippy.core.collection import RegionCollection
from chippy.express.db_query import get_ranked_abs_expr_genes, \
        get_ranked_diff_expr_genes, get_exons
from chippy.util.run_record import RunRecord
from chippy.util.definition import LOG_DEBUG, LOG_INFO, LOG_WARNING, \
    LOG_ERROR, LOG_CRITICAL

__author__ = "Gavin Huttley, Cameron Jack"
__copyright__ = "Copyright 2011-2012, Anuj Pahwa, Gavin Huttley, Cameron Jack"
__credits__ = ["Gavin Huttley", "Cameron Jack"]
__license__ = "GPL"
__maintainer__ = "Cameron Jack"
__email__ = "cameron.jack@anu.edu.au"
__status__ = "Pre-release"
__version__ = '0.1'

class ROI(object):
    def __init__(self, chrom, name, rank, window_start, window_end, strand):
        super(ROI, self).__init__()
        #self.species = species
        self.chrom = chrom
        self.name = name
        self.rank = rank
        self.window_start = window_start # 0-based position
        self.window_end = window_end # 1-based position
        self.strand = strand
        self.counts = numpy.zeros(window_end - window_start)

@display_wrap
def get_counts_ranks_ids(genes, BAMorBED, expr_area,
            max_read_length, count_max_length, window_size=1000,
            rr=RunRecord(), ui=None):
    """ window length=2*window_size (Start of feature is at position 1)
        Build regions of interest (ROI) and return as lists of
        counts, ranks and ensembl_ids, sorted by rank

        max_read_length sets the region that should be counted for any
        read. count_max_length is a boolean for whether to use
        max_read_length or not. CURRENTLY NOT IMPLEMENTED.
    """
    regionsOfInterest = []

    if expr_area.lower() == 'tss':
        for gene in genes:
            start, end, strand = gene.getTssCentredCoords(window_size)
            roi = ROI(gene.chrom, gene.ensembl_id, gene.Rank, start,
                    end, strand)
            regionsOfInterest.append(roi)

    elif expr_area.lower() == 'intron-exon':
        for gene in genes:
            window_list = gene.getAllIntronExonWindows(window_size)
            strand = gene.strand
            for i, window in enumerate(window_list):
                start, end = window
                roi = ROI(gene.chrom, gene.ensembl_id+'_'+str(i),
                        gene.Rank, start, end, strand)
                regionsOfInterest.append(roi)

    if len(regionsOfInterest) == 0:
        rr.display()
        raise RuntimeError('No regions of interest in genome created')

    ROIs, rr = get_region_counts(BAMorBED, regionsOfInterest, rr=rr)
    ROIs = sorted(ROIs, key=lambda roi: roi.rank)
    counts = []; ranks = []; ensembl_ids = []
    for roi in ROIs:
        counts.append(roi.counts)
        ranks.append(roi.rank)
        ensembl_ids.append(roi.name)

    return counts, ranks, ensembl_ids

def centred_counts_for_genes(session, sample_name, expr_area, species,
        BAMorBED, max_read_length, count_max_length, window_size=1000,
        include_target=None, exclude_target=None, rr=RunRecord(),
        test_run=False):
    """returns a RegionCollection object wrapping the counts, ranks etc .."""

    print 'Getting ranked expression instances'
    expressed_genes = get_ranked_abs_expr_genes(session, sample_name,
        include_target=include_target, exclude_target=exclude_target,
        test_run=test_run)

    if not expressed_genes:
        rr.display()
        raise RuntimeError('No expressed genes available for pairing ' +\
                           'with counts data. Halting.')
    rr.addMessage('count_tags', LOG_INFO, 'Sample counts name',
            sample_name)
    rr.addMessage('count_tags', LOG_INFO, 'Total expression data',
         len(expressed_genes))

    print 'Decorating for', len(expressed_genes), 'genes'
    counts, ranks, ensembl_ids = get_counts_ranks_ids(\
        expressed_genes, BAMorBED, expr_area, max_read_length,
        count_max_length, window_size=window_size)

    data = RegionCollection(counts=counts, ranks=ranks,
        labels=ensembl_ids,
        info={'total expressed genes': len(expressed_genes),
                'args': {'window_size': window_size,
                        'max_read_length': max_read_length,
                        'sample_name': sample_name,
                        'species': species}})

    return data, rr

def centred_diff_counts_for_genes(session, sample_name, expr_area, species,
        BAMorBED, max_read_length, count_max_length, window_size,
        multitest_signif_val, include_target=None, exclude_target=None,
        rr=RunRecord(), test_run=False):
    """returns a RegionCollection object wrapping the counts, ranks etc ..
    related to an expression difference experiment"""

    print 'Getting ranked expression difference instances'
    expressed_diff = get_ranked_diff_expr_genes(session, sample_name,
        multitest_signif_val, include_target=include_target,
        exclude_target=exclude_target, test_run=test_run)

    if not expressed_diff:
        rr.display()
        raise RuntimeError('No expressed genes available for pairing ' +\
                'with counts data. Halting.')
    rr.addMessage('count_tags', LOG_INFO, 'Sample diff counts name',
            sample_name)
    rr.addMessage('count_tags', LOG_INFO, 'Total expression data',
            len(expressed_diff))

    counts, ranks, ensembl_ids = get_counts_ranks_ids(\
        expressed_diff, BAMorBED, expr_area, max_read_length,
        count_max_length, window_size=window_size)

    data = RegionCollection(counts=counts, ranks=ranks,
        labels=ensembl_ids,
        info={'total expressed genes': len(expressed_diff),
                'args': {'window_size': window_size,
                        'max_read_length': max_read_length,
                        'sample_name': sample_name,
                        'species': species}})

    return data, rr

