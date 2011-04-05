from sqlalchemy import and_
from sqlalchemy.orm import contains_eager
from sqlalchemy.orm.exc import NoResultFound

from cogent import LoadTable
from cogent.util.progress_display import display_wrap
from cogent.util.misc import flatten

from chippy.express.db_schema import Association, Gene, Transcript, Exon, \
            ExternalGene, Expression, ExpressionDiff, ReferenceFile, Sample, \
            Session, Base, make_session
from chippy.ref.util import chroms
from chippy.util.run_record import RunRecord
from chippy.express.util import single_gene, _one

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

def _get_gene_expression_query(session, ensembl_release, sample_name,
            data_path=None, test=False):
    """returns a query instance"""
    sample = _get_sample(session, sample_name)
    if sample is None:
        raise RuntimeError('Unknown sample name: %s' % sample_name)
    
    if data_path is not None:
        reffile_id = _one(session.query(ReferenceFile.reffile_id).\
                            filter(ReferenceFile.name==data_path))
        if not data_path:
            raise RuntimeError('Unknown data_path %s' % data_path)
        
        reffile_id = reffile_id[0]
    
    if data_path:
        query = session.query(Gene).join(Transcript).filter(
        and_(Transcript.gene_id==Gene.gene_id,
            Gene.ensembl_release==ensembl_release)).join(Association).\
            filter(Transcript.transcript_id==Association.transcript_id).\
            join(Expression).\
            filter(and_(Association.expression_id==Expression.expression_id,
            Expression.sample_id==sample.sample_id,
            Expression.reffile_id==reffile_id)).\
            distinct()
    else:
        query = session.query(Gene).join(Transcript).filter(
        and_(Transcript.gene_id==Gene.gene_id,
            Gene.ensembl_release==ensembl_release)).join(Association).\
            filter(Transcript.transcript_id==Association.transcript_id).\
            join(Expression).\
            filter(and_(Association.expression_id==Expression.expression_id,
            Expression.sample_id==sample.sample_id)).\
            distinct()
    
    return query

def get_samples(session):
    """returns all samples"""
    samples = session.query(Sample).all()
    return samples

def get_transcript_gene_mapping(session, ensembl_release):
    """get ensembl transcript id to chippy gene id mapping"""
    transcript_to_gene = dict(session.query(Transcript.ensembl_id,
        Transcript.gene_id).filter_by(ensembl_release=ensembl_release).all())
    return transcript_to_gene

def get_ensembl_id_transcript_mapping(session, ensembl_release):
    """get ensembl transcript id to chippy transcript id mapping"""
    transcripts = session.query(Transcript.ensembl_id,
                Transcript.transcript_id).filter_by(
                ensembl_release=ensembl_release).all()
    ensembl_id_transcript = dict(transcripts)
    return transcript_to_gene


def get_ranked_expression(session, ensembl_release, sample_name, data_path=None, test=False):
    """returns all ranked genes from a sample"""
    query = _get_gene_expression_query(session, ensembl_release, sample_name,
                    data_path=data_path, test=test)
    genes = query.options(contains_eager('transcripts.expressions')).all()
    genes = sorted((g.getMeanRank(), g) for g in genes)
    return [g for r, g in genes]

def get_ranked_genes_per_chrom(session, ensembl_release, sample_name, chrom, data_path=None, test=False):
    """returns genes from a chromosome"""
    # TODO remove hardcoding for mouse!
    assert chrom in chroms['mouse']
    
    query = _get_gene_expression_query(session, ensembl_release, sample_name,
                             data_path=data_path, test=test)
    
    query = query.filter(Gene.coord_name==chrom)
    genes = query.all()
    genes = sorted((g.getMeanRank(), g) for g in genes)
    return [g for r, g in genes]

def get_external_genes_from_expression_study(session, ensembl_release, external_gene_sample_name, sample_name, test=False):
    external_sample = _get_sample(session, external_gene_sample_name)
    if not external_sample:
        raise RuntimeError('No external_sample with name %s' % \
                        external_gene_sample_name)
    
    query = _get_gene_expression_query(session, ensembl_release, sample_name)
    all = query.all()
    # get the gene id's from the external sample
    external_genes = query.join(ExternalGene).filter(and_(
                    ExternalGene.gene_id==Gene.gene_id,
                    ExternalGene.sample_id==external_sample.sample_id,
                    Gene.ensembl_release==ensembl_release))
    genes = query.all()
    genes = flatten(genes)
    return genes

def get_total_gene_counts(session, ensembl_release, sample_name, data_path=None, test=False):
    """docstring for get_total_gene_counts"""
    query = _get_gene_expression_query(session, ensembl_release, sample_name,
        data_path, test)
    return query.count()


@display_wrap
def diff_expression_study(session, sample_name, data_path, table,
            ensembl_release='58', ensembl_id_label='ENSEMBL', test=False,
            run_record=None, ui=None):
    """returns list of genes and why they're not present in an existing
    expression db"""
    if run_record is None:
        run_record = RunRecord()
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
    
    full_log = LoadTable(header=['Type', 'Id', 'Explanation'], rows=failures,
    title='Difference between db for'\
        +' %(sample)s and %(reffile)s' % dict(reffile=data_path,
                                               sample=sample_name))
    
    # TODO import log constants and change to LOG_INFO etc ..
    failed_probsets = run_record.count("Type == 'probeset'")
    dup_gene_transcripts = run_record.count("Type == 'gene'")
    dberror = run_record.count("Type == 'dberror'")
    run_record.addMessage('diff_expression_study',
        'notice', 'Num genes in db', len(records))
    run_record.addMessage('diff_expression_study',
        'notice', 'Probsets that do not map to a unique gene', failed_probsets)
    run_record.addMessage('diff_expression_study',
        'notice', 'Genes with multiple transcripts', dup_gene_transcripts)
    run_record.addMessage('diff_expression_study',
        'notice', 'Database population errors', dberror)
    
    return run_record, full_log

