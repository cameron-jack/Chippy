"""returns tag counts for a specified collection of genes"""
from __future__ import division

import os

from cogent.util.progress_display import display_wrap

from chippy.core.read_count import WholeChrom
from chippy.core.collection import RegionCollection
from chippy.express.db_query import get_ranked_expression, \
            get_external_genes_from_expression_study
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
def get_count_decorated_expressed_genes(expressed, counts_dir, chrom_names, max_read_length, window_size, ui=None):
    """decorates the Expression instances with a counts attribute, length=2*window_size"""
    # group the expressed by chromosome
    chrom_ordered = {}
    for expressed_gene in ui.series(expressed, noun='Grouping into chromosomes'):
        try:
            chrom_ordered[expressed_gene.gene.coord_name].append(expressed_gene)
        except KeyError:
            chrom_ordered[expressed_gene.gene.coord_name] = [expressed_gene]
    
    assert set(chrom_ordered.keys()) <= set(chrom_names), \
                    'Chromosome mismatch between study and species reference'
    
    n = 0
    total = len(expressed)
    for chrom_name in sorted(chrom_ordered):
        chrom_counts_path = os.path.join(counts_dir,
                    'chr%s.txt.gz' % chrom_name)
        counts = WholeChrom(chrom_counts_path,
                            max_read_length=max_read_length)
        print 'Making full counts array for chromosome %s' % chrom_name
        counts.update()
        
        print '\tDecorating expressed_genes'
        for expressed_gene in chrom_ordered[chrom_name]:
            start, end = expressed_gene.gene.getTssCentredCoords(window_size)
            expressed_gene.counts = counts[start:end].copy()
            n += 1
            if n % 10 == 0:
                ui.display('Decorating genes [%d / %d]' % (n, total), n/total)
        
        del counts
    
    return expressed

@display_wrap
def get_counts_ranks_ensembl_ids(expressed, ui=None):
    """returns separate series for counts, ranks and ensembl_ids"""
    ranks = []
    counts = []
    ensembl_ids = []
    for expressed_gene in ui.series(expressed,
                        noun='Getting counts, ranks and ensembl_ids'):
        ranks.append(expressed_gene.rank)
        counts.append(expressed_gene.counts)
        ensembl_ids.append(expressed_gene.gene.ensembl_id)
    return counts, ranks, ensembl_ids

def _get_decorated_expressed(session, sample_name, species, counts_dir, ensembl_release, max_read_length, window_size, test_run):
    species_chroms = chroms[species]
    
    print 'Getting ranked expression instances'
    expressed = get_ranked_expression(session, ensembl_release, sample_name,
                    test=test_run)
    
    print 'Decorating'
    expressed = get_count_decorated_expressed_genes(expressed, counts_dir,
                            species_chroms, max_read_length, window_size)
    
    for i in range(1, len(expressed)):
        assert expressed[i].rank > expressed[i-1].rank,\
            'Ranks were not sequential'
    
    return expressed

def centred_counts_for_genes(session, sample_name, species, counts_dir,
                    ensembl_release, max_read_length, window_size, test_run):
    """returns a RegionCollection object wrapping the counts, ranks etc .."""
    expressed = _get_decorated_expressed(session, sample_name, species,
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
