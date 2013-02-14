"""returns tag counts for a specified collection of genes"""
from __future__ import division

import sys
sys.path.extend(['..'])
import numpy

from cogent.format.bedgraph import bedgraph

from chippy.core.read_count import get_region_counts
from chippy.core.collection import RegionCollection
from chippy.express.db_query import get_ranked_abs_expr_genes, \
        get_ranked_diff_expr_genes, get_exons
from chippy.util.run_record import RunRecord

__author__ = "Gavin Huttley, Cameron Jack"
__copyright__ = "Copyright 2011-2012, Anuj Pahwa, Gavin Huttley, Cameron Jack"
__credits__ = ["Gavin Huttley", "Cameron Jack"]
__license__ = "GPL"
__maintainer__ = "Cameron Jack"
__email__ = "cameron.jack@anu.edu.au"
__status__ = "Pre-release"
__version__ = '0.1'

class ROI(object):
    def __init__(self, gene, window_start, window_end, label=None):
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
        self.window_start = window_start # 1-based position
        self.window_end = window_end # 1-based position
        try:
            self.counts = numpy.zeros(window_end - window_start, dtype=numpy.uint32)
        except ValueError:
            print self.chrom, self.label, self.rank, self.TSS, \
                    window_start, window_end, self. strand

def write_to_bedgraph(bedgraph_fn, ROIs, rr=RunRecord()):
    """ ROIs sorted by chromosome then start location get written as
    bedgraphs using PyCogent's BEDgraph writer.
    """

    ROIs = sorted(ROIs, key=lambda roi: (roi.chrom, roi.window_start))

    record_tuples = []
    for roi in ROIs:
        for i, c in enumerate(roi.counts):
            if c != 0:
                pos = roi.window_start + i - 1 # bedgraph is 0-offset
                record_tuples.append((roi.chrom, pos, pos, c))
    rr.addInfo('write_to_bedgraph', 'number of records to combine',
            len(record_tuples))
    bedgraph_data = bedgraph(record_tuples, digits=None, name=bedgraph_fn,
            description='Study of filename', color=(0,0,255))

    bgfile = open(bedgraph_fn, 'w')
    bgfile.write(''.join(bedgraph_data)+'\n')
    bgfile.close()

    rr.addInfo('write_to_bedgraph', 'BEDgraph data written to', bedgraph_fn)
    return rr

def get_counts_ranks_ids(genes, BAMorBED, expr_area,
            chr_prefix, window_size=1000, bedgraph=None,
            rr=RunRecord(), ui=None):
    """ window length=2*window_size (Start of feature is at position 1)
        Build regions of interest (ROI) and return as lists of
        counts, ranks and ensembl_ids, sorted by rank
    """
    regionsOfInterest = []

    if expr_area.lower() == 'tss':
        for gene in genes:
            win_start, win_end = gene.getTssCentredCoords(window_size)
            roi = ROI(gene, win_start, win_end)
            regionsOfInterest.append(roi)

    elif expr_area.lower() == 'intron-exon':
        for gene in genes:
            window_list = gene.getAllIntronExonWindows(window_size)
            for i, window in enumerate(window_list):
                win_start, win_end = window
                roi = ROI(gene, win_start, win_end,
                        gene.ensembl_id+'_'+str(i))
                regionsOfInterest.append(roi)

    if not len(regionsOfInterest):
        rr.display()
        raise RuntimeError('No regions of interest in genome created')

    ROIs, rr = get_region_counts(BAMorBED, regionsOfInterest,
            chr_prefix=chr_prefix, rr=rr)

    if bedgraph:
        rr = write_to_bedgraph(bedgraph, ROIs, rr)

    ROIs = sorted(ROIs, key=lambda roi: roi.rank)
    counts = []; ranks = []; ensembl_ids = []
    for i, roi in enumerate(ROIs):
        counts.append(roi.counts)
        ranks.append(roi.rank)
        ensembl_ids.append(roi.gene_id)

    return counts, ranks, ensembl_ids, rr

def centred_counts_for_genes(session, sample_name, expr_area, species,
        BAMorBED, chr_prefix, window_size=1000,
        include_target=None, exclude_target=None,
        bedgraph=None, rr=RunRecord(),
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
    rr.addInfo('count_tags', 'Sample counts name',
            sample_name)
    rr.addInfo('count_tags', 'Total expression data',
         len(expressed_genes))

    print 'Decorating for', len(expressed_genes), 'genes'
    counts, ranks, ensembl_ids, rr = get_counts_ranks_ids(\
            expressed_genes, BAMorBED, expr_area, chr_prefix,
            window_size=window_size, bedgraph=bedgraph, rr=rr)

    data = RegionCollection(counts=counts, ranks=ranks,
        labels=ensembl_ids,
        info={'total expressed genes': len(expressed_genes),
                'args': {'window_size': window_size,
                'sample_name': sample_name,
                'species': species}})

    return data, rr

def centred_diff_counts_for_genes(session, sample_name, expr_area, species,
        BAMorBED, chr_prefix, window_size,
        multitest_signif_val, include_target=None, exclude_target=None,
        bedgraph=None, rr=RunRecord(), test_run=False):
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
    rr.addInfo('count_tags', 'Sample diff counts name',
            sample_name)
    rr.addInfo('count_tags', 'Total expression data',
            len(expressed_diff))

    counts, ranks, ensembl_ids, rr = get_counts_ranks_ids(\
            expressed_diff, BAMorBED, expr_area, chr_prefix,
            window_size=window_size, bedgraph=bedgraph, rr=rr)

    data = RegionCollection(counts=counts, ranks=ranks,
        labels=ensembl_ids,
        info={'total expressed genes': len(expressed_diff),
                'args': {'window_size': window_size,
                'sample_name': sample_name,
                'species': species}})

    return data, rr

