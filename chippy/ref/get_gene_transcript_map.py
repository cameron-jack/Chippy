from cogent import LoadTable
from cogent.db.ensembl import Genome, HostAccount

def gene_transcript_mappings(species, release, account=None, debug=False):
    """returns {gene: [transcripts,,], ..} and {transcript: gene, ..}"""
    transcript_gene = {}
    gene_transcript = {}
    genome = Genome(species, Release=release, account=account)
    bio_types = genome.getDistinct('BioType')
    n = 0
    for bio_type in bio_types:
        genes = genome.getGenesMatching(BioType=bio_type)
        for gene in genes:
            gene_stable_id = gene.StableId
            transcript_ids = [t.StableId for t in gene.Transcripts]
            gene_transcript[gene_stable_id] = transcript_ids
            transcript_gene.update(dict([(gene_stable_id, t)
                        for t in transcript_ids]))
            n += 1
            if n % 100 == 0 and debug:
                return gene_transcript, transcript_gene
    
    return gene_transcript, transcript_gene

if __name__ == "__main__":
    account = HostAccount('cg.anu.edu.au', 'gavin', 'gavin')
    species = 'mouse'
    release = 58
    g, t = gene_transcript_mappings(species, release, account, True)
    
    print g
    print
    print t
