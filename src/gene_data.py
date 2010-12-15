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
        self.order = None
        self.stable_id = stable_id
    
    def __str__(self):
        return '%s(chrom=%s, index=%s)' % (self.stable_id, self.chrom, self.index)
    
    def __repr__(self):
        return '%s(%s, %s)' % (self.stable_id, self.chrom, self.index)
    
    def __cmp__(self, other):
        return cmp((self.chrom, self.index), (other.chrom, other.index))

class GetGeneIndexes(object):
    """class to store  / return gene indices"""
    def __init__(self):
        super(GetGeneIndexes, self).__init__()
        gene_indices = LoadTable(os.path.join(data_dir,
                        'mouse_gene_coords.txt'), sep='\t')
        names_indices = []
        for stable_id, chrom_name, index in gene_indices.getRawData(['StableId',
                            'CoordName', 'Index']):
            gene = GeneIndex(stable_id, chrom_name, index)
            names_indices += [(stable_id, gene)]
        
        self._gene_table = gene_indices
        self.names_indices = dict(names_indices)
    
    def _call(self, gene_id):
        return self.names_indices[gene_id]
    
    __call__ = _call
    
    @property
    def GeneTable(self):
        return self._gene_table
    

GeneIndexes = GetGeneIndexes()

def get_upstream_coords(infile_name, window_size):
    """Given a file with gene information (EnsembleID, start, end, strand),
    return a data structure with EnsembleID, chromosome, TSS, coordinate for
    window_size upstream from TSS, and the strand that the gene is found on."""

    table = LoadTable(infile_name, sep='\t')
    table = table.sorted(columns=['CoordName', 'Start'])
    rows = []
    for row in table:
        strand = row['Strand']
        if strand == -1:
            TSS = row['End']
            start_coord = TSS
            end_coord = TSS + window_size
        else:
            TSS = row['Start']
            start_coord = TSS - window_size
            end_coord = TSS

        rows += [[row['StableId'], row['CoordName'], start_coord, end_coord,
                  row['Strand']]]
    return rows

if __name__ == "__main__":

    from light_seq import LightSeq
    from cogent.db.ensembl import HostAccount, Species, Genome

    # get the coordinates of interest
    coords = get_upstream_coords('../data/mouse_gene_coords.txt', 500)

    # Ensembl connect information
    account = HostAccount('cg.anu.edu.au','compgen','compgen')
    specie = Genome(Species='Mouse', Release = 58, account=account)

    outfile_name = '../data/mouse_TATA.fasta'
    outfile = open(outfile_name, 'w')

    for [stable_id, chrom, start, end, strand] in coords:
        print 'Processing %s on Chromosome %s'%(stable_id, str(chrom))
        gene = specie.getRegion(CoordName=str(chrom), Start=start, End=end)
        if strand == -1:
            seq_str = str(gene.Seq.rc())
        else:
            seq_str = str(gene.Seq)

        seq = LightSeq(Seq=seq_str, Name=stable_id)
        outfile.write(seq.toFasta() + '\n')

    outfile.close()