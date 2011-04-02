from sqlalchemy import and_
from sqlalchemy.orm import contains_eager
from sqlalchemy.orm.exc import NoResultFound

from cogent.util.progress_display import display_wrap
from cogent.util.misc import flatten

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

def get_samples(session):
    """returns all samples"""
    samples = session.query(Sample).all()
    return samples

def get_ranked_expression(session, ensembl_release, sample_name, test=False):
    """returns all ranked genes from a sample"""
    sample = _get_sample(session, sample_name)
    if sample is None:
        raise RuntimeError('Unknown sample name: %s' % sample_name)
    
    expressed = session.query(Expression).join(Gene).filter(
            and_(Gene.ensembl_release==ensembl_release,
                 Expression.sample_id==sample.sample_id)).order_by(
                 Expression.rank).options(contains_eager('gene')).all()
    return expressed

def get_ranked_genes_per_chrom(session, ensembl_release, sample_name, chrom, test=False):
    """returns genes from a chromosome"""
    # TODO remove hardcoding for mouse!
    assert chrom in chroms['mouse']
    sample = _get_sample(session, sample_name)
    if sample is None:
        raise RuntimeError('Unknown sample name: %s' % sample_name)
    
    expressed = session.query(Expression).join(Gene).filter(
            and_(Gene.coord_name==chrom,
                 Gene.ensembl_release==ensembl_release,
                 Expression.sample_id==sample.sample_id)).order_by(
                 Expression.rank).options(contains_eager('gene')).all()
    
    return expressed

def get_external_genes_from_expression_study(session, ensembl_release, external_gene_sample_name, sample_name, test=False):
    expression_sample = _get_sample(session, sample_name)
    external_sample = _get_sample(session, external_gene_sample_name)
    if expression_sample is None or external_sample is None:
        raise RuntimeError('Unknown sample name(s): %s / %s' % (
                                external_gene_sample_name, sample_name))
    
    # get the gene id's from the external sample
    genes = session.query(Gene.gene_id).join(ExternalGene).filter(and_(
                    ExternalGene.gene_id==Gene.gene_id,
                    ExternalGene.sample_id==external_sample.sample_id,
                    Gene.ensembl_release==ensembl_release)).all()
    genes = flatten(genes)
    if not genes:
        raise RuntimeError('No external genes found from %s' % external_sample)
    
    # get the expressed instances from joining to the expression sample
    expressed = session.query(Expression).join(Gene).filter(
            and_(Expression.gene_id.in_(genes),
                 Expression.sample_id==expression_sample.sample_id)).order_by(
                 Expression.rank).options(contains_eager('gene')).all()
    return expressed

def get_total_gene_counts(session, ensembl_release, sample_name, test=False):
    """docstring for get_total_gene_counts"""
    sample = _get_sample(session, sample_name)
    if sample is None:
        raise RuntimeError('Unknown sample name: %s' % sample_name)
    
    total = session.query(Expression).join(Gene).filter(
            and_(Gene.ensembl_release==ensembl_release,
                 Expression.sample_id==sample.sample_id)).count()
    return total
