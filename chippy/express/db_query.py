from sqlalchemy import and_
from sqlalchemy.orm import contains_eager
from sqlalchemy.orm.exc import NoResultFound

from cogent import LoadTable
from cogent.util.progress_display import display_wrap
from cogent.util.misc import flatten

from chippy.express.db_schema import Gene, Transcript, Exon, \
            ExternalGene, Expression, ExpressionDiff, ReferenceFile, Sample, \
            Session, Base, make_session
from chippy.express.db_populate import get_transcript_gene_mapping, single_gene
from chippy.ref.util import chroms

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"
__version__ = '0.1'

def _one(query):
    """returns result if found, False otherwise"""
    try:
        result = query.one()
    except NoResultFound:
        result = False
    
    return result

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

def get_ranked_expression(session, ensembl_release, sample_name, data_path=None, test=False):
    """returns all ranked genes from a sample"""
    sample = _get_sample(session, sample_name)
    if data_path is not None:
        reffile = _one(session.query(ReferenceFile).filter_by(name=data_path))
        if reffile is None:
            raise RuntimeError('No reference file record for that path')
    
    if sample is None:
        raise RuntimeError('Unknown sample name: %s' % sample_name)
    if test:
        if data_path is None:
            total = session.query(Expression).filter(
            Expression.sample_id==sample.sample_id).count()
        else:
            total = session.query(Expression).filter(
            and_(Expression.sample_id==sample.sample_id,
                Expression.reference_file_id==reffile.reference_file_id)
                ).count()
        
        print 'Num expression records for sample in db', total
        
    if data_path is None:
        query = session.query(Expression).join(Gene).filter(
            and_(Gene.ensembl_release==ensembl_release,
                 Expression.sample_id==sample.sample_id)).order_by(
                 Expression.rank).options(contains_eager('gene'))
    else:
        query = session.query(Expression).join(Gene).filter(
            and_(Gene.ensembl_release==ensembl_release,
                 Expression.sample_id==sample.sample_id,
                 Expression.reference_file_id==reffile.reference_file_id)
                 ).order_by(Expression.rank).options(contains_eager('gene'))
    
    expressed = query.all()
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


@display_wrap
def diff_expression_study(session, sample_name, data_path, table,
            ensembl_release='58', ensembl_id_label='ENSEMBL', test=False, ui=None):
    """returns list of genes and why they're not present in an existing
    expression db"""
    # query db for every gene linked to this expression study
    records = get_ranked_expression(session, ensembl_release, sample_name,
            data_path=data_path, test=test)
    # 
    transcript_to_gene = get_transcript_gene_mapping(session, ensembl_release)
    failures = []
    kept_gene_ids_probesets = []
    for record in ui.series(table, noun='Adding expression diffs'):
        transcript_ids = record[ensembl_id_label]
        gene_id = single_gene(transcript_to_gene, transcript_ids)
        probeset_id = record['probeset_id']
        if gene_id is None:
            failures.append(['probeset', probeset_id, 
                            'probeset id maps to multiple genes'])
            continue
        kept_gene_ids_probesets.append([gene_id, probeset_id])
    
    in_db_gene_ids = set([e.gene_id for e in records])
    in_table_kept_gene_ids = set([g for g,p in kept_gene_ids_probesets])
    in_db_not_table = in_db_gene_ids - in_table_kept_gene_ids
    if in_db_not_table:
        failures.extend([('dberror', gid,
                    'population error! Contact developers')
                for gid in in_db_not_table])
    
    in_table_not_db = in_table_kept_gene_ids - in_db_gene_ids
    if in_table_not_db:
        failures.extend([('gene', gid, 'gene mapped to multiple probesets')
                                                for gid in in_table_not_db])
    
    run_record = LoadTable(header=['Type', 'Id', 'Explanation'], rows=failures,
    title='Difference between db for'\
        +' %(sample)s and %(reffile)s' % dict(reffile=data_path,
                                               sample=sample_name))
    
    failed_probsets = run_record.count("Type == 'probeset'")
    dup_gene_transcripts = run_record.count("Type == 'gene'")
    dberror = run_record.count("Type == 'dberror'")
    summary = LoadTable(header=['Label', 'Number'],
        rows=[['Num genes in db', len(records)],
              ['Probsets that do not map to a unique gene', failed_probsets],
              ['Genes with multiple transcripts', dup_gene_transcripts],
              ['Database population errors', dberror]],
              title='Summary of difference')
    
    return summary, run_record

