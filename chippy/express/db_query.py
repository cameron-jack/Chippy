from sqlalchemy import and_
from sqlalchemy.orm.exc import NoResultFound

from cogent.util.progress_display import display_wrap

from chippy.express.db_schema import Gene, Transcript, Exon, \
            ExternalGene, Expression, ExpressionDiff, ReferenceFile, Sample, \
            Session, Base, make_session
from chippy.ref.util import chroms

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"
__version__ = '0.1'

def _get_sample(session, sample_name):
    try:
        sample = session.query(Sample).filter(Sample.name==sample_name).one()
    except NoResultFound:
        sample = None
    return sample

def get_ranked_expression(session, ensembl_release, sample_name, test=False):
    """returns all ranked genes from a sample"""
    sample = _get_sample(session, sample_name)
    if sample is None:
        raise RuntimeError('Unknown sample name: %s' % sample_name)
    
    expressed = session.query(Expression).join(Gene).filter(
            and_(Gene.ensembl_release==ensembl_release,
                 Expression.sample_id==sample.sample_id)).order_by(
                 Expression.rank).all()
    return expressed

def get_ranked_genes_per_chrom(session, ensembl_release, sample_name, chrom, test=False):
    """returns genes from a chromosome"""
    assert chrom in chroms['mouse']
    sample = _get_sample(session, sample_name)
    if sample is None:
        raise RuntimeError('Unknown sample name: %s' % sample_name)
    
    expressed = session.query(Expression).join(Gene).filter(
            and_(Gene.coord_name==chrom,
                 Gene.ensembl_release==ensembl_release,
                 Expression.sample_id==sample.sample_id)).order_by(
                 Expression.rank).all()
    
    return expressed
