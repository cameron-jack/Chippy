"""given an rdump, find the within chromosome gene indexes for the listed
genes"""
from __future__ import division
from cogent import LoadTable
from parse_r_dump import RDumpToTable

# we first get the rdumped genes and the ensembl transcript to gene id table
def get_gene_ids(rdump_filename):
    """returns gene IDs ordered by their rank in the rdump file"""
    rdump = RDumpToTable(rdump_filename)
    rdump = rdump.filtered(lambda x: '---' not in x, columns='ENSEMBL')
    transcript_ids = rdump.getRawData('ENSEMBL')
    transcript_to_gene = LoadTable('../data/ensembl_transcript_gene_ids.txt',sep='\t')
    transcript_to_gene = transcript_to_gene.getRawData(['TranscriptId', 'GeneId'])
    transcript_to_gene = dict(transcript_to_gene)
    rows = []
    problem = 0
    for index, transcript_group in enumerate(transcript_ids):
        gene_ids = set()
        for transcript_id in transcript_group:
            try:
                gene_id = transcript_to_gene[transcript_id]
            except KeyError:
                continue
            gene_ids.update([gene_id])
        
        if len(gene_ids) != 1:
            problem += 1
            print '\nThe following transcript ids: %s\nMatch multiple genes: %s' % (str(transcript_group), gene_ids)
            continue
        rows += [list(gene_ids)[0]]
    print 'There were %d problem transcripts:' % (problem)
    return rows

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

def get_gene_indices(rdump_filename):
    gene_ids = get_gene_ids(rdump_filename)
    get_gene_index = GetGeneIndexes()
    genes = [get_gene_index(gene_id) for gene_id in gene_ids]
    data = sorted([(gene.chrom, gene) for gene in genes])
    return [g for c,g in data]


if __name__ == "__main__":
    gene_ids = get_gene_ids('../tests/data/rdump_sample.txt')
    get_gene_index = GetGeneIndexes()
    genes = [get_gene_index(gene_id) for gene_id in gene_ids]
    print genes
    