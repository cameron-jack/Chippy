from cogent import LoadTable
from cogent.db.ensembl import *

account = HostAccount('cg.anu.edu.au', 'compgen', 'compgen')
mouse = Genome('mouse', Release=58,account=account)

def run():
    transcript_gene = []
    bio_types = mouse.getDistinct('BioType')
    for bio_type in bio_types:
        genes = mouse.getGenesMatching(BioType=bio_type)
        for index, gene in enumerate(genes):
            print index, gene.Symbol
            gene_stable_id = gene.StableId
            gene_symbol = gene.Symbol
            for transcript in gene.Transcripts:
                transcript_stable_id = transcript.StableId
                transcript_gene += [[bio_type, transcript_stable_id, gene_stable_id, gene_symbol]]


    table = LoadTable(header=['BioType', 'TranscriptId', 'GeneId', 'GeneSymbol'], rows=transcript_gene)
    table.writeToFile('../data/ensembl_transcript_gene_ids.txt', sep='\t')

if __name__ == "__main__":
    run()
