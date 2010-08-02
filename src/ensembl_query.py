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
    chr_MT = Mouse.getRegion(CoordName = 19)
    out = open('chr_19.fasta','w')
    out.write(chr_19.Seq.toFasta()+'\n')
    out.close()
    
ensembl_query()









