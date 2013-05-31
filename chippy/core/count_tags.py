"""returns tag counts for a specified collection of genes"""
from __future__ import division

import sys
sys.path.extend(['..'])

from cogent.format.bedgraph import bedgraph
from chippy.core.region_of_interest import ROI
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

def write_to_bedgraph(bedgraph_fn, ROIs):
    """ ROIs sorted by chromosome then start location get written as
        bedgraphs using PyCogent's BEDgraph writer.
    """
    rr = RunRecord('write_to_bedgraph')
    ROIs = sorted(ROIs, key=lambda roi: (roi.chrom, roi.start))

    record_tuples = []
    for roi in ROIs:
        for i, c in enumerate(roi.counts):
            if c != 0:
                pos = roi.start + i # BEDgraph is 0-based like Python
                record_tuples.append((roi.chrom, pos, pos, c))
    rr.addInfo('number of records to combine', len(record_tuples))
    bedgraph_data = bedgraph(record_tuples, digits=None, name=bedgraph_fn,
            description='Study of filename', color=(0,0,255))

    bgfile = open(bedgraph_fn, 'w')
    bgfile.write(''.join(bedgraph_data)+'\n')
    bgfile.close()

    rr.addInfo('BEDgraph data written to', bedgraph_fn)

def get_counts_ranks_ids(genes, BAMorBED, expr_area,
            chr_prefix, window_radius=1000, bedgraph=None, ui=None):
    """ window length=2*window_radius (Start of feature is right of centre.
        Build regions of interest (ROI) and return as lists of
        counts, ranks and ensembl_ids, sorted by rank.
        All coordinates in Python 0-offset space.

        Returns counts, ranks, labels, and possible normalisation factors:
        number of tags (reads), number of bases counted.
    """
    rr = RunRecord('get_counts_ranks_ids')
    regionsOfInterest = []

    if expr_area.lower() == 'tss':
        for gene in genes:
            win_start, win_end = gene.getTssCentredCoords(window_radius)
            roi = ROI(gene, win_start, win_end)
            regionsOfInterest.append(roi)

    elif expr_area.lower() == 'intron-exon':
        for gene in genes:
            window_list = gene.getAllIntronExonWindows(window_radius)
            for i, window in enumerate(window_list):
                win_start, win_end = window
                roi = ROI(gene, win_start, win_end,
                        gene.ensembl_id+'_'+str(i))
                regionsOfInterest.append(roi)

    if not len(regionsOfInterest):
        rr.display()
        raise RuntimeError('No regions of interest in genome created')

    ROIs, num_tags, num_bases = get_region_counts(BAMorBED,
            regionsOfInterest, chr_prefix=chr_prefix)

    if bedgraph:
        write_to_bedgraph(bedgraph, ROIs)

    ROIs = sorted(ROIs, key=lambda roi: roi.rank)
    counts = []; ranks = []; ensembl_ids = []
    for i, roi in enumerate(ROIs):
        counts.append(roi.counts)
        ranks.append(roi.rank)
        ensembl_ids.append(roi.gene_id)

    return counts, ranks, ensembl_ids, num_tags, num_bases

def centred_counts_for_genes(session, sample_name, expr_area, species,
        BAMorBED, chr_prefix, window_radius=1000,
        include_target=None, exclude_target=None,
        bedgraph=None, test_run=False):
    """returns a RegionCollection object wrapping the counts, ranks etc .."""
    rr = RunRecord('centred_counts_for_genes')

    print 'Getting ranked expression instances'
    expressed_genes = get_ranked_abs_expr_genes(session, sample_name,
            include_target=include_target, exclude_target=exclude_target,
            test_run=test_run)

    if not expressed_genes:
        rr.display()
        raise RuntimeError('No expressed genes available for pairing ' +\
                           'with counts data. Halting.')
    rr.addInfo('Sample counts name', sample_name)
    rr.addInfo('Total expression data', len(expressed_genes))

    print 'Decorating for', len(expressed_genes), 'genes'
    counts, ranks, ensembl_ids, num_tags, num_bases =\
            get_counts_ranks_ids(expressed_genes, BAMorBED, expr_area,
            chr_prefix, window_radius=window_radius, bedgraph=bedgraph)

    data = RegionCollection(counts=counts, ranks=ranks,
            labels=ensembl_ids,
            info={'total expressed genes': len(expressed_genes),
                'args': {'window_radius': window_radius,
                'sample_name': sample_name,
                'species': species,
                'tag count': num_tags,
                'base count': num_bases}})

    return data

def centred_diff_counts_for_genes(session, sample_name, expr_area, species,
        BAMorBED, chr_prefix, window_radius,
        multitest_signif_val, include_target=None, exclude_target=None,
        bedgraph=None, test_run=False):
    """returns a RegionCollection object wrapping the counts, ranks etc ..
    related to an expression difference experiment"""
    rr = RunRecord('centred_diff_counts_for_genes')

    print 'Getting ranked expression difference instances'
    expressed_diff = get_ranked_diff_expr_genes(session, sample_name,
            multitest_signif_val, include_target=include_target,
            exclude_target=exclude_target, test_run=test_run)

    if not expressed_diff:
        rr.display()
        raise RuntimeError('No expressed genes available for pairing ' +\
                'with counts data. Halting.')
    rr.addInfo('Sample diff counts name', sample_name)
    rr.addInfo('Total expression data', len(expressed_diff))

    counts, ranks, ensembl_ids, num_tags, num_bases =\
            get_counts_ranks_ids(expressed_diff, BAMorBED, expr_area,
            chr_prefix, window_radius=window_radius, bedgraph=bedgraph)

    data = RegionCollection(counts=counts, ranks=ranks,
        labels=ensembl_ids,
        info={'total expressed genes': len(expressed_diff),
                'args': {'window_radius': window_radius,
                'sample_name': sample_name,
                'species': species,
                'tag count': num_tags,
                'base count': num_bases}})

    return data

