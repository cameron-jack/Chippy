"""returns tag counts for a specified collection of genes"""
from __future__ import division

import sys
sys.path.extend(['..'])
import numpy

from cogent.util.progress_display import display_wrap

from chippy.core.read_count import get_combined_counts, read_all_beds
from chippy.core.collection import RegionCollection
from chippy.express.db_query import get_ranked_abs_expr_genes, \
        get_ranked_diff_expr_genes, get_ranked_genes_per_chrom, \
        get_diff_ranked_genes_per_chrom, get_chroms
from chippy.util.util import grouped_by_chrom
from chippy.util.run_record import RunRecord
from chippy.util.definition import LOG_DEBUG, LOG_INFO, LOG_WARNING, \
    LOG_ERROR, LOG_CRITICAL

__author__ = "Gavin Huttley, Cameron Jack"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley, Cameron Jack"
__credits__ = ["Gavin Huttley, Cameron Jack"]
__license__ = "GPL"
__maintainer__ = "Cameron Jack"
__email__ = "cameron.jack@anu.edu.au"
__status__ = "Pre-release"
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

class ROI(object):
    def __init__(self, chrom, name, window_start, window_end, strand):
        super(ROI, self).__init__()
        #self.species = species
        self.chrom = chrom
        self.name = name
        self.window_start = window_start
        self.window_end = window_end
        self.strand = strand
        self.count_array = numpy.zeros(window_end - window_start)

@display_wrap
def get_count_decorated_expressed_genes(genes, BAMorBED, expr_area,
            max_read_length, count_max_length, window_size=1000,
            rr=RunRecord(), ui=None):
    """ decorates the Expression instances with a counts attribute,
        length=2*window_size (Start of feature is at position 1)

    Get genes -> get locations, get windows, for whichever location type
    call read_counts(BAMorBED, RegionsOfInterest)

    """
    regionsOfInterest = []

    if expr_area.lower() == 'tss':
        for gene in genes:
            start, end, strand = gene.getTssCentredCoords(window_size)
            roi = ROI(gene.chrom, gene.ensembl_id, start, end, strand)
            regionsOfInterest.append(ROI)
    elif expr_area.lower() == 'intron-exon':
        for gene in genes:
            intron_window_list = gene.getAllIntron3primeWindows(window_size)
            strand = gene.strand



            # strand is being used here as a stride so if minus strand
            # then start > end and with stride==-1 we reverse the counts so
            # they are all 5' to 3'
            feature_counts = counts[start:end:strand].copy()
            feature = Feature(feature_counts, gene.Rank, gene.ensembl_id)
            features.append(feature)

            if intron_window_list is not None:
                for intron_window in intron_window_list:
                    start = intron_window[0]
                    end = intron_window[1]
                    feature_counts = counts[start:end:strand].copy()
                    feature = Feature(feature_counts, gene.Rank, gene.ensembl_id)
                    features.append(feature)

    n = 0
    total = len(genes)
    summed_counts = {}
    features = [] # this is the output set

    counts = get_combined_counts(BAMorBED,
            max_read_length, count_max_length)
    if counts == None:
            rr.display()
            raise RuntimeError('No data loaded from ' + BAMorBED)
    summed_counts[chrom_name] = counts.counts.sum()
    print '\tGetting read counts for genes'
    # All below needs to be above to first get the regions required and create
    # ROIs, which can then be sent to the BAM or BED readers
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
                raise RuntimeError('Count tags, expression area invalid: ' + expr_area)

            n += 1
            if n % 10 == 0:
                ui.display('Getting counts for gene [' + str(n) + ', ' + \
                        str(total) + ' / ' + str(n/total) + '%]')
        del counts
    
    return features, summed_counts

def _get_decorated_expressed(session, sample_name, expr_area, species,
        BAMorBED, max_read_length, count_max_length, window_size,
        include_target=None, exclude_target=None, rr=RunRecord(),
        test_run=False):

    print 'Getting ranked expression instances'
    expressed_genes = get_ranked_abs_expr_genes(session, sample_name,
            include_target=include_target, exclude_target=exclude_target,
            test_run=test_run)

    if not expressed_genes:
        rr.display()
        raise RuntimeError('No expressed genes available for pairing ' +\
                'with counts data. Halting.')

    print 'Decorating for %d genes' % len(expressed)
    expressed, summed_counts = get_count_decorated_expressed_genes(\
            expressed_genes, BAMorBED, expr_area, max_read_length,
            count_max_length, window_size=window_size)
    
    return expressed, summed_counts

def _get_decorated_expressed_diff(session, sample_name, expr_area, species,
        BAMorBED, max_read_length, count_max_length, window_size,
        multitest_signif_val, include_target=None, exclude_target=None,
        test_run=False):
    """ Combine counts and expression data """
    
    print 'Getting ranked expression difference instances'
    expressed_diff = get_ranked_diff_expr_genes(session, sample_name,
            multitest_signif_val, include_target=include_target,
            exclude_target=exclude_target, test_run=test_run)

    if not expressed_diff:
        rr.display()
        raise RuntimeError('No expressed genes available for pairing ' +\
                           'with counts data. Halting.')

    print 'Decorating'
    expressed_diff, summed_counts = get_count_decorated_expressed_genes(\
            expressed_diff, BAMorBED, expr_area, max_read_length,
            count_max_length, window_size=window_size)

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

def centred_counts_for_genes(session, sample_name, expr_area, species,
        BAMorBED, max_read_length, count_max_length, window_size=1000,
        include_target=None, exclude_target=None, run_record=None,
        test_run=False):
    """returns a RegionCollection object wrapping the counts, ranks etc .."""
    
    expressed, summed_counts = _get_decorated_expressed(session, sample_name,
            expr_area, species, BAMorBED, max_read_length,
            count_max_length, window_size=window_size,
            include_target=include_target, exclude_target=exclude_target,
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
        BAMorBED, max_read_length, count_max_length, window_size,
        multitest_signif_val, include_target=None, exclude_target=None,
        run_record=None, test_run=False):
    """returns a RegionCollection object wrapping the counts, ranks etc ..
    related to an expression difference experiment"""

    if run_record is None:
        run_record = RunRecord()

    expressed_diff, summed_counts = _get_decorated_expressed_diff(session,
            sample_name, expr_area, species, BAMorBED, max_read_length,
            count_max_length, window_size, multitest_signif_val, include_target,
            exclude_target, test_run)
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

