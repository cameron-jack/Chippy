"""basic class to store gene coordinates plus a factory function to return
these instances for a collection of Ensembl stable id's"""
import os

from cogent import LoadTable

from util import data_dir

class GeneIndex(object):
    """basic object to store a chromosome and an index"""
    def __init__(self, stable_id, chrom, index):
        super(GeneIndex, self).__init__()
        self.chrom = chrom
        self.index = index
        self.stable_id = stable_id
    
    def __str__(self):
        return '%s(chrom=%s, index=%s)' % (self.stable_id, self.chrom, self.index)
    
    def __repr__(self):
        return '%s(%s, %s)' % (self.stable_id, self.chrom, self.index)

def GetGeneIndexes():
    """returns func to get indexes for groups of genes"""
    gene_indices = LoadTable(os.path.join(data_dir,
                    'mouse_gene_coords.txt'), sep='\t')
    names_indices = []
    for stable_id, chrom_name, index in gene_indices.getRawData(['StableId',
                        'CoordName', 'Index']):
        gene = GeneIndex(stable_id, chrom_name, index)
        names_indices += [(stable_id, gene)]
    
    names_indices = dict(names_indices)
    def call(gene_id):
        return names_indices[gene_id]
    
    return call
