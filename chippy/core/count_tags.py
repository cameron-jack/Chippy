"""returns tag counts for a specified collection of genes"""
from __future__ import division

import sys
sys.path.extend(['..'])

from cogent.format.bedgraph import bedgraph
from chippy.core.region_of_interest import ROI
from chippy.core.read_count import get_region_counts
from chippy.core.collection import RegionCollection
from chippy.express.db_query import get_genes_by_ranked_expr, \
        get_genes_by_ranked_diff, get_species
from chippy.util.run_record import RunRecord
from chippy.express.util import sample_types

__author__ = 'Cameron Jack, Gavin Huttley'
__copyright__ = 'Copyright 2011-2013, Gavin Huttley, Cameron Jack, Anuj Pahwa'
__credits__ = ['Gavin Huttley', 'Cameron Jack']
__license__ = 'GPL'
__maintainer__ = 'Cameron Jack'
__email__ = 'cameron.jack@anu.edu.au'
__status__ = 'Pre-release'
__version__ = '0.3'

def write_to_bedgraph(bedgraph_fn, ROIs):
    """
        ROIs sorted by chromosome then start location get written as
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

def get_counts_ranks_ids(genes, BAMorBED, feature_type,
            chr_prefix, window_upstream, window_downstream,
            bedgraph=None, ui=None):
    """
        Build regions of interest (ROI) and return as lists of counts, ranks
        and ensembl_ids, sorted by rank.

        All coordinates in Python 0-offset space.

        Returns counts, ranks, labels, and possible normalisation factors:
        number of tags (reads), number of bases counted.
    """
    rr = RunRecord('get_counts_ranks_ids')
    regionsOfInterest = []

    if feature_type.lower() == 'tss':
        # Transcriptional Start Tite
        for gene in genes:
            win_start, win_end = gene.getTssWindowCoords(\
                    window_upstream, window_downstream)
            roi = ROI(gene, win_start, win_end)
            regionsOfInterest.append(roi)

    elif feature_type.lower() == 'utr_exon':
        # Boundary of the 5' UTR with the first Exon
        # Except when overlapping with 5' TSS
        for gene in genes:
            win_start, win_end = gene.getUTRExonWindowCoords(\
                    window_upstream, window_downstream)
            if win_start is not None and win_end is not None:
                roi = ROI(gene, win_start, win_end)
                regionsOfInterest.append(roi)


    elif feature_type.lower() == 'intron_exon':
        # All intron/exon boundaries, except when overlapping UTR/Exon boundary
        for gene in genes:
            window_list = gene.getIntronExonWindowCoords(\
                    window_upstream, window_downstream)
            for i, window in enumerate(window_list):
                win_start, win_end = window
                if win_start is not None and win_end is not None:
                    roi = ROI(gene, win_start, win_end,
                            gene.ensembl_id+'_'+str(i))
                    regionsOfInterest.append(roi)

    elif feature_type.lower() == 'exon_intron':
        # all exon/intron boundaries, except when overlapping UTR/Exon boundary
        for gene in genes:
            window_list = gene.getExonIntronWindowCoords(\
                    window_upstream, window_downstream)
            for i, window in enumerate(window_list):
                win_start, win_end = window
                if win_start is not None and win_end is not None:
                    roi = ROI(gene, win_start, win_end,
                            gene.ensembl_id+'_'+str(i))
                    regionsOfInterest.append(roi)

    elif feature_type.lower() == 'exon_utr':
        # Boundary of the last exon with the 3' UTR - except when overlapping
        # with 3' end of gene.
        for gene in genes:
            win_start, win_end = gene.getExonUTRWindowCoords(\
                    window_upstream, window_downstream)
            if win_start is not None and win_end is not None:
                roi = ROI(gene, win_start, win_end)
                regionsOfInterest.append(roi)

    elif feature_type.lower() == 'gene_3p':
        # 3' end of gene
        for gene in genes:
            win_start, win_end = gene.getTssWindowCoords(\
                    window_upstream, window_downstream)
            roi = ROI(gene, win_start, win_end)
            regionsOfInterest.append(roi)

    else:
        rr.dieOnCritical('Feature type not supported', feature_type)

    if not len(regionsOfInterest):
        rr.dieOnCritical('No Regions of Interest created')

    ROIs, num_tags, num_bases, mapped_tags = get_region_counts(BAMorBED,
            regionsOfInterest, chr_prefix=chr_prefix)

    if bedgraph:
        write_to_bedgraph(bedgraph, ROIs)

    ROIs = sorted(ROIs, key=lambda roi: roi.rank)
    counts = []; ranks = []; ensembl_ids = []
    for i, roi in enumerate(ROIs):
        counts.append(roi.counts)
        ranks.append(roi.rank)
        ensembl_ids.append(roi.gene_id)

    return counts, ranks, ensembl_ids, num_tags, num_bases, mapped_tags

def counts_for_genes(session, sample_name, sample_type,
        feature_type, BAMorBED, chr_prefix, window_upstream,
        window_downstream, include_target=None, exclude_target=None,
        bedgraph=None, multitest_signif_val=None):
    """returns a RegionCollection object wrapping the counts, ranks etc .."""
    rr = RunRecord('counts_for_genes')

    expressed_genes = None
    if sample_type == sample_types['exp_absolute']:
        print 'Getting ranked expression instances'
        expressed_genes = get_genes_by_ranked_expr(session, sample_name,
                include_target=include_target, exclude_target=exclude_target)

    elif sample_type == sample_types['exp_diff']:
        print 'Getting ranked expression difference instances'
        expressed_genes = get_genes_by_ranked_diff(session, sample_name,
                multitest_signif_val=multitest_signif_val,
                include_target=include_target, exclude_target=exclude_target)

    else:
        rr.dieOnCritical('Sample type not supported', sample_type)

    if expressed_genes is None:
        rr.dieOnCritical('Expressed genes', 'not present')

    rr.addInfo('Sample counts name', sample_name)
    rr.addInfo('Total expression data', len(expressed_genes))

    print 'Getting counts of region', feature_type,'for', \
            len(expressed_genes), 'genes'
    counts, ranks, ensembl_ids, num_tags, num_bases, mapped_tags =\
            get_counts_ranks_ids(expressed_genes, BAMorBED, feature_type,
            chr_prefix, window_upstream, window_downstream,
            bedgraph=bedgraph)

    data = RegionCollection(counts=counts, ranks=ranks,
            labels=ensembl_ids,
            info={'total expressed genes': len(expressed_genes),
                'args': {
                    'window_upstream': window_upstream,
                    'window_downstream': window_downstream,
                    'feature_type': feature_type,
                    'sample_name': sample_name,
                    'species': get_species(session),
                    'tag count': num_tags,
                    'base count': num_bases,
                    'mapped tags': mapped_tags
                }
            })

    return data
