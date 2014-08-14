from sqlalchemy.orm.exc import NoResultFound

__author__ = "Gavin Huttley, Cameron Jack"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__maintainer__ = "Cameron Jack"
__email__ = "cameron.jack@anu.edu.au"
__status__ = "Pre-release"
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

sample_types = {
    'abs_expr' : 'Absolute expression data',
    'diff_expr' : 'Differential expression data',
    'target_genes' : 'Target gene list'
}

sample_type_names = [key for key in sorted(sample_types.keys())]
sample_type_desc = [val for val in sorted(sample_types.values())]
