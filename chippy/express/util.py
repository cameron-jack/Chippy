__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"
__version__ = '0.1'


def single_gene(transcript_to_gene, transcript_ids):
    gene_ids = set()
    for tid in transcript_ids:
        try:
            gid = transcript_to_gene[tid]
        except KeyError:
            continue
        gene_ids.update([gid])
    if len(gene_ids) == 1:
        value = gid
    else:
        value = None
    return value

def _one(query):
    """returns result if found, False otherwise"""
    try:
        result = query.one()
    except NoResultFound:
        result = False
    
    return result
