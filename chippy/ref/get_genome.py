from cogent.db.ensembl import Genome, HostAccount

def get_chrom_seqs(species, release, chroms, account=None):
    """yields sequence objects for the indicated chromosomes from Ensembl"""
    genome = Genome(species, Release=release, account=account)
    for chrom in chroms:
        region = genome.getRegion(CoordName = chrom)
        seq = region.Seq
        name = 'chr_%s.fasta' % chrom
        seq.Name = name
        yield seq


if __name__ == "__main__":
    host = 'cg.anu.edu.au'
    username = password = 'gavin'
    release = 58
    chroms = map(str, range(1,20)+['X', 'Y'])
    account = HostAccount(host, username, password)
    print account
    species = 'mouse'
    for chrom in get_chrom_seqs(species, release, chroms, account):
        print chrom.Name, len(chrom)
        fasta = chrom[:10].toFasta()
        print fasta
        # outfile_name = os.path.join(out_dir, )
        # outfile = open(outfile_name, 'w')
        # outfile.write(fasta+'\n')
        # outfile.close()
