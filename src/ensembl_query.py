import os
from cogent import LoadTable
from cogent.db.ensembl import HostAccount, Species, Genome

try:
    username, password = os.environ['ENSEMBL_ACCOUNT'].split()
except KeyError:
    print 'Need ENSEMBL_ACCOUNT environment variable, quitting!'
    exit(-1)

release = 58
chroms = map(str, range(1,20)+['X', 'Y'])
account = HostAccount('cg.anu.edu.au', username, password)
mouse = Genome(Species='Mouse', Release = release, account=account)

def sample_ensembl(out_dir, write_seqs=False):
    for chrom in chroms:
        region = mouse.getRegion(CoordName = chrom)
        if write_seqs:
            outfile_name = os.path.join(out_dir, 'chr_%s.fasta' % chrom)
            outfile = open(outfile_name, 'w')
            out.write(region.Seq.toFasta()+'\n')
            out.close()
    

def get_chrom_lengths(out_dir):
    """writes a table of chromosome lengths"""
    rows = []
    for chrom in chroms:
        region = mouse.getRegion(CoordName = chrom)
        rows += [[chrom, len(region)]]
    table = LoadTable(header=['chrom', 'length'], rows=rows)
    print table
    table.writeToFile(os.path.join(out_dir,
                    'mouse_chrom_lengths_release_%d.txt' % release), sep='\t')

if __name__ == "__main__":
    # sample_ensembl('', release=58, write_seqs=False)
    get_chrom_lengths("../data/")