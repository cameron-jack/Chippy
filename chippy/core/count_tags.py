"""returns tag counts for a specified collection of genes"""
from __future__ import division

import os, sys

from cogent.util.progress_display import display_wrap

from chippy.core.read_count import get_combined_counts, read_all_beds
from chippy.core.collection import RegionCollection
from chippy.express.db_query import get_ranked_expression, \
        get_ranked_expression_diff, get_ranked_genes_per_chrom, \
        get_diff_ranked_genes_per_chrom
from chippy.ref.util import chroms
from chippy.util.util import grouped_by_chrom
from chippy.util.run_record import RunRecord
from chippy.util.definition import LOG_DEBUG, LOG_INFO, LOG_WARNING, \
    LOG_ERROR, LOG_CRITICAL

__author__ = "Gavin Huttley, Cameron Jack"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley, Cameron Jack"
__credits__ = ["Gavin Huttley, Cameron Jack"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"
__version__ = '0.1'

class Feature:
    """ Abstraction for gene, exon, intron and any other genomic feature """
    def __init__(self, counts=None, ranks=None, ids=None):
        self.counts = counts
        self.Rank = ranks
        self.ensembl_id = ids

    # This is required to support sorting
    def __repr__(self):
        return repr((self.counts, self.Rank, self.ensembl_id))

@display_wrap
def get_count_decorated_expressed_genes(genes, counts_dir, expr_area, chrom_names,
            max_read_length, count_max_length, window_size=1000, ui=None):
    """decorates the Expression instances with a counts attribute, length=2*window_size"""

    # group the genes by chromosome
    chrom_ordered = grouped_by_chrom(genes)
    
    assert set(chrom_ordered.keys()) <= set(chrom_names), \
                    'Chromosome mismatch between study and species reference'

    # import all bed data first and then mine it per chrom later
    bed_reps = read_all_beds(counts_dir)

    n = 0
    total = len(genes)
    summed_counts = {}
    features = [] # this is the output set
    chrom_names = [c for l, c in sorted([(len(c), c) for c in chrom_ordered])]

    for chrom_name in chrom_names:
        print 'Making full counts array for chromosome %s' % chrom_name
        counts = get_combined_counts(counts_dir, bed_reps, chrom_name, max_read_length,
                                    count_max_length)
        if counts == None:
            continue
        summed_counts[chrom_name] = counts.counts.sum()
        print '\tGetting read counts for genes'
        for gene in chrom_ordered[chrom_name]:
            if expr_area.lower() == 'tss':
                start, end, strand = gene.getTssCentredCoords(window_size)
                # strand is being used here as a stride so if minus strand
                # then start > end and with stride==-1 we reverse the counts so
                # they are all 5' to 3'
                feature_counts = counts[start:end:strand].copy()
                feature = Feature(feature_counts, gene.Rank, gene.ensembl_id)
                features.append(feature)

            elif expr_area.lower() == 'intron_3p':
                # intron 3' and exon 5' regions - except first exon 5'
                intron_window_list = gene.getAllIntron3primeWindows(window_size)
                strand = gene.strand
                if intron_window_list is not None:
                    for intron_window in intron_window_list:
                        start = intron_window[0]
                        end = intron_window[1]
                        feature_counts = counts[start:end:strand].copy()
                        feature = Feature(feature_counts, gene.Rank, gene.ensembl_id)
                        features.append(feature)

            elif expr_area.lower() == 'exon_3p':
                # exon 3' and intron 5', also includes final exon 3'
                exon_window_list = gene.getAllExon3primeWindows(window_size)
                strand = gene.strand
                if exon_window_list is not None:
                    for exon_window in exon_window_list:
                        start = exon_window[0]
                        end = exon_window[1]
                        feature_counts = counts[start:end:strand].copy()
                        feature = Feature(feature_counts, gene.Rank, gene.ensembl_id)
                        features.append(feature)

            elif expr_area.lower() == 'both_3p':
                # Both intron and exon 3' regions, also includes final exon 3'
                strand = gene.strand
                intron_window_list = gene.getAllIntron3primeWindows(window_size)
                if intron_window_list is not None:
                    for intron_window in intron_window_list:
                        start = intron_window[0]
                        end = intron_window[1]
                        feature_counts = counts[start:end:strand].copy()
                        feature = Feature(feature_counts, gene.Rank, gene.ensembl_id)
                        features.append(feature)
                exon_window_list = gene.getAllExon3primeWindows(window_size)
                if exon_window_list is not None:
                    for exon_window in exon_window_list:
                        start = exon_window[0]
                        end = exon_window[1]
                        feature_counts = counts[start:end:strand].copy()
                        feature = Feature(feature_counts, gene.Rank, gene.ensembl_id)
                        features.append(feature)

            else:
                raise RuntimeError('Count tags, expression area invalid: %s' % (expr_area))

            n += 1
            if n % 10 == 0:
                ui.display('Getting counts for gene [%d / %d]' % \
                                                    (n, total), n/total)
        del counts
    
    return features, summed_counts

def _get_decorated_expressed(session, sample_name, expr_area, species, chrom,
        counts_dir, max_read_length, count_max_length, window_size,
        include_genes=None, exclude_genes=None, test_run=False):
    species_chroms = chroms[species]

    msg = 'Getting ranked expression instances%s'
    if chrom is None:
        print msg % ''
        expressed = get_ranked_expression(session, sample_name,
                include_genes=include_genes, exclude_genes=exclude_genes,
                test_run=test_run)
    else:
        print msg % (' for chrom ' + chrom)
        expressed = get_ranked_genes_per_chrom(session,
                sample_name, chrom, include_genes=include_genes,
                exclude_genes=exclude_genes, test_run=test_run)
    
    print 'Decorating for %d genes' % len(expressed)
    expressed, summed_counts = get_count_decorated_expressed_genes(expressed,
            counts_dir, expr_area, species_chroms, max_read_length, count_max_length,
            window_size=window_size)
    
    return expressed, summed_counts

def _get_decorated_expressed_diff(session, sample_name, expr_area, species,
        chrom, counts_dir, max_read_length, count_max_length, window_size,
        multitest_signif_val, include_genes=None, exclude_genes=None,
        test_run=False):
    
    species_chroms = chroms[species]

    msg = 'Getting ranked expression difference instances%s'
    if chrom is None:
        print msg % ''
        expressed_diff = get_ranked_expression_diff(session, sample_name,
                multitest_signif_val, include_genes=include_genes,
                exclude_genes=exclude_genes, test_run=test_run)
    else:
        print msg % (' for chrom ' + chrom)
        expressed_diff = get_diff_ranked_genes_per_chrom(session, sample_name,
                multitest_signif_val, chrom, include_genes=include_genes,
                exclude_genes=exclude_genes, test_run=test_run)

    print 'Decorating'
    expressed_diff, summed_counts = get_count_decorated_expressed_genes(expressed_diff,
            counts_dir, expr_area, species_chroms, max_read_length, count_max_length,
            window_size=window_size)

    return expressed_diff, summed_counts

@display_wrap
def get_counts_ranks_ensembl_ids(features, ui=None):
    """returns separate series for counts, ranks and ensembl_ids"""
    ranks = []
    counts = []
    ensembl_ids = []
    # This sort is required for the Feature class to match the Gene class
    features = sorted(features, key=lambda features: features.Rank)
    for feature in ui.series(features, noun='Getting counts, ranks and ensembl_ids'):
        ranks.append(feature.Rank)
        counts.append(feature.counts)
        ensembl_ids.append(feature.ensembl_id)
    return counts, ranks, ensembl_ids

def centred_counts_for_genes(session, sample_name, expr_area, species, chrom,
        counts_dir, max_read_length, count_max_length, window_size=1000,
        include_genes=None, exclude_genes=None, run_record=None,
        test_run=False):
    """returns a RegionCollection object wrapping the counts, ranks etc .."""
    
    expressed, summed_counts = _get_decorated_expressed(session, sample_name,
            expr_area, species, chrom, counts_dir, max_read_length,
            count_max_length, window_size=window_size,
            include_genes=include_genes, exclude_genes=exclude_genes,
            test_run=test_run)
    total_expressed_genes = len(expressed)
    run_record.addMessage('count_tags', LOG_INFO, 'Sample counts name',
                          sample_name)
    run_record.addMessage('count_tags', LOG_INFO, 'Total expression data',
                          total_expressed_genes)

    if total_expressed_genes == 0:
        return None, run_record

    counts, ranks, ensembl_ids = get_counts_ranks_ensembl_ids(expressed)

    data = RegionCollection(counts=counts, ranks=ranks,
        labels=ensembl_ids,
        info={'total expressed genes': total_expressed_genes,
                'args': {'window_size': window_size,
                        'max_read_length': max_read_length,
                        'sample_name': sample_name,
                        'species': species,
                        'summed_counts': summed_counts}})

    return data, run_record

def centred_diff_counts_for_genes(session, sample_name, expr_area, species,
        chrom, counts_dir, max_read_length, count_max_length, window_size,
        multitest_signif_val, include_genes=None, exclude_genes=None,
        run_record=None, test_run=False):
    """returns a RegionCollection object wrapping the counts, ranks etc ..
    related to an expression difference experiment"""

    if run_record is None:
        run_record = RunRecord()

    expressed_diff, summed_counts = _get_decorated_expressed_diff(session,
            sample_name, expr_area, species, chrom, counts_dir, max_read_length,
            count_max_length, window_size, multitest_signif_val, include_genes,
            exclude_genes, test_run)
    total_expressed_diff_genes = len(expressed_diff)
    
    run_record.addMessage('count_tags', LOG_INFO, 'Sample diff counts name',
                          sample_name)
    run_record.addMessage('count_tags', LOG_INFO, 'Total expression data',
                          total_expressed_diff_genes)

    if total_expressed_diff_genes == 0:
        return None, run_record

    counts, ranks, ensembl_ids = get_counts_ranks_ensembl_ids(expressed_diff)

    data = RegionCollection(counts=counts, ranks=ranks,
        labels=ensembl_ids,
        info={'total expressed genes': total_expressed_diff_genes,
                'args': {'window_size': window_size,
                        'max_read_length': max_read_length,
                        'sample_name': sample_name,
                        'species': species,
                        'summed_counts': summed_counts}})

    return data, run_record

def centred_counts_target_genes(session, target_genes_sample_name, 
    sample_name, species, counts_dir, max_read_length,
    window_size):
    """returns RegionCollection object for target gene list sample from
    the sample_name"""
    expressed = get_target_genes_from_expression_study(session,
                target_genes_sample_name, sample_name)
    chrom_names = chroms[species]
    expressed, summed_counts = get_count_decorated_expressed_genes(expressed,
        counts_dir, chrom_names, max_read_length, window_size)
    total_expressed_genes = len(expressed)
    
    counts, ranks, ensembl_ids = get_counts_ranks_ensembl_ids(expressed)
    
    data = RegionCollection(counts=counts, ranks=ranks,
        labels=ensembl_ids,
        info={'total expressed genes': total_expressed_genes,
                'args': {'window_size': window_size,
                        'max_read_length': max_read_length,
                        'sample_name': sample_name,
                        'species': species,
                        'summed_counts': summed_counts}})
    
    return data

