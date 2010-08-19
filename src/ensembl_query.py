import os
from cogent import LoadTable
from cogent.db.ensembl import HostAccount, Species, Genome

try:
    username, password = os.environ['ENSEMBL_ACCOUNT'].split()
except KeyError:
    print 'Need ENSEMBL_ACCOUNT environment variable, quitting!'
    exit(-1)

def sample_ensembl(out_dir, release, write_seqs=False):
    account = HostAccount('cg.anu.edu.au', username, password)
    
    mouse = Genome(Species='Mouse', Release = release, account=account)
    
    chroms = map(str, range(1,20)+['X', 'Y'])
    for chrom in chroms:
        region = mouse.getRegion(CoordName = chrom)
        chrom_lengths += [[chrom, len(region)]]
        
        if write_seqs:
            outfile_name = os.path.join(out_dir, 'chr_%s.fasta' % chrom)
            outfile = open(outfile_name, 'w')
            out.write(region.Seq.toFasta()+'\n')
            out.close()
    

if __name__ == "__main__":
    sample_ensembl('', release=58, write_seqs=False)