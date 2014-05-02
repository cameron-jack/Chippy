"""returns tag counts for a specified collection of genes"""
from __future__ import division

import sys
sys.path.extend(['..'])

from cogent.format.bedgraph import bedgraph
from chippy.core.region_of_interest import ROI
from chippy.core.read_count import get_region_counts
from chippy.core.collection import RegionCollection
from chippy.express import db_query
from chippy.util.run_record import RunRecord
from chippy.express.util import sample_types

from chippy.util.definition import NULL_STRAND, PLUS_STRAND, MINUS_STRAND

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
        current_count = 0
        for i, c in enumerate(roi.counts):
            if c != 0:
                if current_count != c and current_count == 0:
                    current_count = c
                elif current_count != c:
                    if roi.strand == PLUS_STRAND:
                        record_tuples.append((roi.chrom, roi.start, pos, c))
                    else:
                        record_tuples.append((roi.chrom, roi.end, pos, c))

                if roi.strand == PLUS_STRAND:
                    pos = roi.start + i # BEDgraph is 0-based like Python
                else:
                    pos = roi.end + i # but is unidirectional
        del roi

    rr.addInfo('number of records to combine', len(record_tuples))
    bedgraph_data = bedgraph(record_tuples, digits=None, name=bedgraph_fn,
            description='Study of filename', color=(0,0,255))

    bgfile = open(bedgraph_fn, 'w')
    bgfile.write(''.join(bedgraph_data)+'\n')
    bgfile.close()

    rr.addInfo('BEDgraph data written to', bedgraph_fn)

def write_BED_windows(BED_windows_fn, regions_of_interest):
    """
        Create a single BED3 file entry for each ROI
    """
    rr = RunRecord('write_BED_windows')
    rr.addInfo('number of records to write', len(regions_of_interest))

    ROIs = sorted(regions_of_interest, key=lambda x:\
            (x.chrom, x.start, x.end))

    with open(BED_windows_fn, 'w') as bed:
        for roi in ROIs:
            chrom = roi.chrom
            start = roi.start
            end = roi.end # no BED indexing adjustment needed
            line = '\t'.join(map(str, [chrom, start, end])) + '\n'
            bed.write(line)

    rr.addInfo('BED windows written to', BED_windows_fn)

def get_counts_ranks_ids(genes, BAMorBED, feature_type, chr_prefix,
        window_upstream, window_downstream, bedgraph_fn=None,
        BED_windows_fn=None, chrom_size=300000000, ui=None):
    """
        Build regions of interest (ROI) and return as lists of counts, ranks
        and ensembl_ids, sorted by rank.

        All coordinates in Python 0-offset space.

        Returns counts, ranks, labels, and possible normalisation factors:
        number of tags (reads), number of bases counted.
    """
    rr = RunRecord('get_counts_ranks_ids')
    regions_of_interest = []

    if feature_type.lower() == 'tss':
        # Transcriptional Start Tite
        for gene in genes:
            win_start, win_end = gene.getTssWindowCoords(\
                    window_upstream, window_downstream)
            roi = ROI(gene, win_start, win_end)
            regions_of_interest.append(roi)

    elif feature_type.lower() == 'utr_exon':
        # Boundary of the 5' UTR with the first Exon
        # Except when overlapping with 5' TSS
        for gene in genes:
            win_start, win_end = gene.getUTRExonWindowCoords(\
                    window_upstream, window_downstream)
            if win_start is not None and win_end is not None:
                roi = ROI(gene, win_start, win_end)
                regions_of_interest.append(roi)

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
                    regions_of_interest.append(roi)

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
                    regions_of_interest.append(roi)

    elif feature_type.lower() == 'exon_utr':
        # Boundary of the last exon with the 3' UTR - except when overlapping
        # with 3' end of gene.
        for gene in genes:
            win_start, win_end = gene.getExonUTRWindowCoords(\
                    window_upstream, window_downstream)
            if win_start is not None and win_end is not None:
                roi = ROI(gene, win_start, win_end)
                regions_of_interest.append(roi)

    elif feature_type.lower() == 'gene_3p':
        # 3' end of gene
        for gene in genes:
            win_start, win_end = gene.getGene3PrimeWindowCoords(\
                    window_upstream, window_downstream)
            roi = ROI(gene, win_start, win_end)
            regions_of_interest.append(roi)

    else:
        rr.dieOnCritical('Feature type not supported', feature_type)

    if not len(regions_of_interest):
        rr.dieOnCritical('No Regions of Interest created')

    if BED_windows_fn:
        write_BED_windows(BED_windows_fn, regions_of_interest)

    ROIs, num_tags, num_bases, mapped_tags = get_region_counts(BAMorBED,
            regions_of_interest, chr_prefix=chr_prefix, chrom_size=chrom_size)

    if bedgraph_fn:
        write_to_bedgraph(bedgraph_fn, ROIs)

    ROIs = sorted(ROIs, key=lambda roi: roi.rank)
    counts = []; ranks = []; ensembl_ids = []
    for i, roi in enumerate(ROIs):
        counts.append(roi.counts)
        ranks.append(roi.rank)
        ensembl_ids.append(roi.gene_id)

    return counts, ranks, ensembl_ids, num_tags, num_bases, mapped_tags

def counts_for_genes(session, sample_name, feature_type, BAMorBED, chr_prefix,
        window_upstream, window_downstream, include_target=None,
        exclude_target=None, bedgraph_fn=None, multitest_signif_val=None,
        BED_windows_fn=None, chrom_size=300000000):
    """returns a RegionCollection object wrapping the counts, ranks etc .."""
    rr = RunRecord('counts_for_genes')

    expressed_genes = None
    if sample_name in db_query.get_expr_sample_names(session):
        print 'Getting ranked expression instances'
        expressed_genes = db_query.get_genes_by_ranked_expr(session, sample_name,
            include_target=include_target, exclude_target=exclude_target)
    elif sample_name in db_query.get_diff_sample_names(session):
        print 'Getting ranked expression difference instances'
        expressed_genes = db_query.get_genes_by_ranked_diff(session, sample_name,
            multitest_signif_val=multitest_signif_val,
            include_target=include_target, exclude_target=exclude_target)
    else:
        rr.dieOnCritical('Sample must be either', 'Absolute or Differential')

    if expressed_genes is None:
        rr.dieOnCritical('Expressed genes', 'not present')

    rr.addInfo('Sample counts name', sample_name)
    rr.addInfo('Total expression data', len(expressed_genes))

    print 'Getting counts of region', feature_type,'for', \
            len(expressed_genes), 'genes'
    counts, ranks, ensembl_ids, num_tags, num_bases, mapped_tags =\
            get_counts_ranks_ids(expressed_genes, BAMorBED, feature_type,
            chr_prefix, window_upstream, window_downstream,
            bedgraph_fn=bedgraph_fn, BED_windows_fn=BED_windows_fn,
            chrom_size=chrom_size)

    data = RegionCollection(counts=counts, ranks=ranks,
            labels=ensembl_ids,
            info={'total expressed genes': len(expressed_genes),
                'args': {
                    'window_upstream': window_upstream,
                    'window_downstream': window_downstream,
                    'feature_type': feature_type,
                    'sample_name': sample_name,
                    'species': db_query.get_species(session),
                    'tag count': num_tags,
                    'base count': num_bases,
                    'mapped tags': mapped_tags
                }
            })

    return data
