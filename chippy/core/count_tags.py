"""returns tag counts for a specified collection of genes"""
from __future__ import division

import os

from cogent.util.progress_display import display_wrap

from chippy.core.read_count import WholeChrom
from chippy.core.collection import RegionCollection
from chippy.express.db_query import get_ranked_expression, \
        get_ranked_genes_per_chrom, get_external_genes_from_expression_study
from chippy.ref.util import chroms

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"
__version__ = '0.1'

@display_wrap
def get_count_decorated_expressed_genes(genes, counts_dir, chrom_names, max_read_length, window_size, ui=None):
    """decorates the Expression instances with a counts attribute, length=2*window_size"""
    # group the genes by chromosome
    chrom_ordered = {}
    for gene in ui.series(genes, noun='Grouping into chromosomes'):
        try:
            chrom_ordered[gene.coord_name].append(gene)
        except KeyError:
            chrom_ordered[gene.coord_name] = [gene]
    
    assert set(chrom_ordered.keys()) <= set(chrom_names), \
                    'Chromosome mismatch between study and species reference'
    
    n = 0
    total = len(genes)
    for chrom_name in sorted(chrom_ordered):
        chrom_counts_path = os.path.join(counts_dir,
                    'chr%s.txt.gz' % chrom_name)
        counts = WholeChrom(chrom_counts_path,
                            max_read_length=max_read_length)
        print 'Making full counts array for chromosome %s' % chrom_name
        counts.update()
        
        print '\tDecorating genes'
        for gene in chrom_ordered[chrom_name]:
            start, end, strand = gene.getTssCentredCoords(window_size)
            gene.counts = counts[start:end:strand].copy()
            n += 1
            if n % 10 == 0:
                ui.display('Decorating genes [%d / %d]' % (n, total), n/total)
        
        del counts
    
    return genes

@display_wrap
def get_counts_ranks_ensembl_ids(genes, ui=None):
    """returns separate series for counts, ranks and ensembl_ids"""
    ranks = []
    counts = []
    ensembl_ids = []
    for gene in ui.series(genes, noun='Getting counts, ranks and ensembl_ids'):
        ranks.append(gene.MeanScore)
        counts.append(gene.counts)
        ensembl_ids.append(gene.ensembl_id)
    return counts, ranks, ensembl_ids

def _get_decorated_expressed(session, sample_name, species, chrom, counts_dir, ensembl_release, max_read_length, window_size, test_run):
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
    expressed = get_count_decorated_expressed_genes(expressed, counts_dir,
                            species_chroms, max_read_length, window_size)
    
    return expressed

def centred_counts_for_genes(session, sample_name, species, chrom, counts_dir,
                    ensembl_release, max_read_length, window_size, test_run):
    """returns a RegionCollection object wrapping the counts, ranks etc .."""
    expressed = _get_decorated_expressed(session, sample_name, species, chrom,
        counts_dir, ensembl_release, max_read_length, window_size, test_run)
    total_expressed_genes = len(expressed)
    
    counts, ranks, ensembl_ids = get_counts_ranks_ensembl_ids(expressed)
    
    data = RegionCollection(counts=counts, ranks=ranks,
        labels=ensembl_ids,
        info={'total expressed genes': total_expressed_genes,
                'args': {'window_size': window_size,
                        'max_read_length': max_read_length,
                        'ensembl_release': ensembl_release,
                        'sample_name': sample_name,
                        'species': species}})
    
    return data

def centred_counts_external_genes(session, external_genes_sample_name, 
    sample_name, species, counts_dir, ensembl_release, max_read_length,
    window_size):
    """returns RegionCollection object for external gene list sample from
    the sample_name"""
    expressed = get_external_genes_from_expression_study(session,
                ensembl_release, external_genes_sample_name, sample_name)
    chrom_names = chroms[species]
    expressed = get_count_decorated_expressed_genes(expressed, counts_dir,
        chrom_names, max_read_length, window_size)
    total_expressed_genes = len(expressed)
    
    counts, ranks, ensembl_ids = get_counts_ranks_ensembl_ids(expressed)
    
    data = RegionCollection(counts=counts, ranks=ranks,
        labels=ensembl_ids,
        info={'total expressed genes': total_expressed_genes,
                'args': {'window_size': window_size,
                        'max_read_length': max_read_length,
                        'ensembl_release': ensembl_release,
                        'sample_name': sample_name,
                        'species': species}})
    
    return data

if __name__ == "__main__":
    pass
