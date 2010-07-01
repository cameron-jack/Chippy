import os
import sqlalchemy as sql
from cogent.db.ensembl import HostAccount, Species, Genome
def ensembl_query():
    
    if 'ENSEMBL_ACCOUNT' in os.environ:
        username, password = os.environ['ENSEMBL_ACCOUNT'].split()
        account = HostAccount('cg.anu.edu.au','compgen','compgen')
    else:
        account = None
    Mouse = Genome(Species='Mouse', Release = 58, account=account)
    chr_MT = Mouse.getRegion(CoordName = chrMT)
    out = open('chr_MT.fasta','w')
    out.write(chr_MT.Seq.toFasta()+'\n')
    out.close()
    
ensembl_query()









