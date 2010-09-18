"""basic class to store gene coordinates plus a factory function to return
these instances for a collection of Ensembl stable id's"""
from cogent import LoadTable

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
    gene_indices = LoadTable('../data/mouse_gene_coords.txt', sep='\t')
    
    names_indices = []
    for chrom_name in range(1,20)+['X', 'Y']:
        chrom = gene_indices.filtered(lambda x: x == chrom_name, columns='CoordName')
        names = chrom.getRawData('StableId')
        for i in range(len(names)):
            gene = GeneIndex(names[i], chrom_name, i)
            names_indices += [(names[i], gene)]
    
    names_indices = dict(names_indices)
    def call(gene_id):
        return names_indices[gene_id]
    
    return call
