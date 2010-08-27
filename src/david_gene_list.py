import os
import sqlalchemy as sql
from cogent.db.ensembl import HostAccount, Species, Genome

account = HostAccount('cg.anu.edu.au','compgen','compgen')

mouse = Genome(Species='Mouse', Release = 58, account=account)
genes = mouse.getGenesMatching(Symbol='E2f1')
lines = []
for gene in genes:
    stab_id = gene.StableId
    symbol = gene.Symbol
    data = '\t'.join([symbol, stab_id])
fh = open('test_list_gene.txt', 'a')
fh.write(data + '\n')
