"""returns tag counts for a specified collection of genes"""
from __future__ import division

import os, sys

from cogent.util.progress_display import display_wrap

from chippy.core.read_count import WholeChrom, get_combined_counts
from chippy.core.collection import RegionCollection
from chippy.express.db_query import get_ranked_expression, \
        get_ranked_genes_per_chrom
from chippy.ref.util import chroms
from chippy.util.util import grouped_by_chrom

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"
__version__ = '0.1'

@display_wrap
def get_count_decorated_expressed_genes(genes, counts_dir, chrom_names, max_read_length, count_max_length, window_size, ui=None):
    """decorates the Expression instances with a counts attribute, length=2*window_size"""
    # group the genes by chromosome
    chrom_ordered = grouped_by_chrom(genes)
    
    assert set(chrom_ordered.keys()) <= set(chrom_names), \
                    'Chromosome mismatch between study and species reference'
    
    n = 0
    total = len(genes)
    summed_counts = {}
    chrom_names = [c for l, c in sorted([(len(c), c) for c in chrom_ordered])]
    for chrom_name in chrom_names:
        print 'Making full counts array for chromosome %s' % chrom_name
        counts = get_combined_counts(counts_dir, chrom_name, max_read_length,
                                    count_max_length)
        summed_counts[chrom_name] = counts.counts.sum()
        print '\tGetting read counts for genes'
        for gene in chrom_ordered[chrom_name]:
            start, end, strand = gene.getTssCentredCoords(window_size)
            # strand is being used here as a stride so if minus strand
            # then start > end and with stride==-1 we reverse the counts so
            # they are all 5' to 3'
            gene.counts = counts[start:end:strand].copy()
            n += 1
            if n % 10 == 0:
                ui.display('Getting counts for gene [%d / %d]' % \
                                                    (n, total), n/total)
        
        del counts
    
    return genes, summed_counts

@display_wrap
def get_counts_ranks_ensembl_ids(genes, ui=None):
    """returns separate series for counts, ranks and ensembl_ids"""
    ranks = []
    counts = []
    ensembl_ids = []
    for gene in ui.series(genes, noun='Getting counts, ranks and ensembl_ids'):
        ranks.append(gene.Rank)
        counts.append(gene.counts)
        ensembl_ids.append(gene.ensembl_id)
    return counts, ranks, ensembl_ids

def _get_decorated_expressed(session, sample_name, species, chrom, counts_dir, ensembl_release, max_read_length, count_max_length, window_size, test_run):
    species_chroms = chroms[species]
    
    msg = 'Getting ranked expression instances%s'
    if chrom is None:
        print msg % ''
        expressed = get_ranked_expression(session, ensembl_release,
                    sample_name, test_run=test_run)
    else:
        print msg % (' for chrom ' + chrom)
        expressed = get_ranked_genes_per_chrom(session, ensembl_release,
                    sample_name, chrom, test_run=test_run)
    
    print 'Decorating'
    expressed, summed_counts = get_count_decorated_expressed_genes(expressed,
            counts_dir, species_chroms, max_read_length, count_max_length,
            window_size)
    
    return expressed, summed_counts

def centred_counts_for_genes(session, sample_name, species, chrom, counts_dir,
                    ensembl_release, max_read_length, count_max_length,
                    window_size, test_run):
    """returns a RegionCollection object wrapping the counts, ranks etc .."""
    expressed, summed_counts = _get_decorated_expressed(session, sample_name,
        species, chrom, counts_dir, ensembl_release, max_read_length,
        count_max_length, window_size, test_run)
    total_expressed_genes = len(expressed)
    if total_expressed_genes == 0:
        sys.stderr.write("No expression data for '%s'\n" % sample_name)
        return
    
    counts, ranks, ensembl_ids = get_counts_ranks_ensembl_ids(expressed)
    
    data = RegionCollection(counts=counts, ranks=ranks,
        labels=ensembl_ids,
        info={'total expressed genes': total_expressed_genes,
                'args': {'window_size': window_size,
                        'max_read_length': max_read_length,
                        'ensembl_release': ensembl_release,
                        'sample_name': sample_name,
                        'species': species,
                        'summed_counts': summed_counts}})
    
    return data

def centred_counts_external_genes(session, external_genes_sample_name, 
    sample_name, species, counts_dir, ensembl_release, max_read_length,
    window_size):
    """returns RegionCollection object for external gene list sample from
    the sample_name"""
    expressed = get_external_genes_from_expression_study(session,
                ensembl_release, external_genes_sample_name, sample_name)
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
                        'ensembl_release': ensembl_release,
                        'sample_name': sample_name,
                        'species': species,
                        'summed_counts': summed_counts}})
    
    return data

