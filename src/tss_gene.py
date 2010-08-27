import os
import sqlalchemy as sql
from cogent.db.ensembl import HostAccount, Species, Genome

account = HostAccount('cg.anu.edu.au','compgen','compgen')

mouse = Genome(Species='Mouse', Release = 58, account=account)
genes = mouse.getGenesMatching(BioType='protein_coding')
count = 0
lines = []
for gene in genes:
    stab_id = gene.StableId
    gene_loc = str(gene.Location.CoordName)
    start = str(gene.Location.Start)
    end = str(gene.Location.End)
    strand = str(gene.Location.Strand)
    lines.append ('\t'.join([stab_id, gene_loc, start, end, strand]))
    
    count += 1
    if count == 100:
        break
data = '\n'.join(lines)
recs = open('gene_records', 'w')
recs.write(data)
recs.close()

