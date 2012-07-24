from sqlalchemy import and_
from sqlalchemy.orm import contains_eager, joinedload
from sqlalchemy.orm.exc import NoResultFound
import pickle
import gzip
from cogent import LoadTable
from cogent.util.progress_display import display_wrap
from cogent.util.misc import flatten

from chippy.express.db_schema import Gene, Exon, \
            TargetGene, Expression, ExpressionDiff, ReferenceFile, Sample, \
            Session, Base, make_session
from chippy.ref.util import chroms
from chippy.util.run_record import RunRecord
from chippy.express.util import single_gene, _one

__author__ = "Gavin Huttley, Cameron Jack"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley, Cameron Jack"
__credits__ = ["Gavin Huttley, Cameron Jack"]
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

def get_gene_expression_query(session, sample_name, data_path=None,
        test_run=False):
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
        # used to reconstruct the origin of a sample
        query = session.query(Expression).join(Gene).filter(
                and_(Expression.sample_id==sample.sample_id,
                Expression.reffile_id==reffile_id)).\
                options(contains_eager('gene'))
    else:
        query = session.query(Expression).join(Gene).filter(
                Expression.sample_id==sample.sample_id).\
                options(contains_eager('gene'))
    
    return query

def get_gene_expression_diff_query(session, sample_name, data_path=None,
        multitest_signif_val=None, test_run=False):
    """returns a query instance for expression difference """
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
        if multitest_signif_val is None:
            query = session.query(ExpressionDiff).join(Gene).filter(
                and_(ExpressionDiff.sample_id==sample.sample_id,
                    ExpressionDiff.reffile_id==reffile_id)).\
                    options(contains_eager('gene'))
        else:
            query = session.query(ExpressionDiff).join(Gene).filter(
                and_(ExpressionDiff.sample_id==sample.sample_id,
                    ExpressionDiff.multitest_signif==multitest_signif_val,
                    ExpressionDiff.reffile_id==reffile_id)).\
                    options(contains_eager('gene'))
    else:
        if multitest_signif_val is None:
            query = session.query(ExpressionDiff).join(Gene).filter(
                and_(ExpressionDiff.sample_id==sample.sample_id)).\
                    options(contains_eager('gene'))
        else:
            query = session.query(ExpressionDiff).join(Gene).filter(
                and_(ExpressionDiff.sample_id==sample.sample_id,
                    ExpressionDiff.multitest_signif==multitest_signif_val)).\
                    options(contains_eager('gene'))

    return query

def get_samples(session):
    """returns all samples"""
    samples = session.query(Sample).all()
    return samples

def get_target_sample(session):
    """returns samples for target genes"""
    query = session.query(Sample).join(TargetGene).distinct()
    return query.all()

def get_genes(session, chrom=None, biotype='protein_coding', stable_ids=None):
    """returns the Gene's for the indicated release, chrom, biotype or ensembl stable ID.
    
    Note: if stable_ids provided, all arguments other are ignored.
    """
    if type(stable_ids) == str:
        stable_ids = [stable_ids]
    
    if stable_ids:
        condition = Gene.ensembl_id.in_(stable_ids)
    elif chrom and biotype:
        condition = and_(Gene.coord_name==str(chrom), Gene.biotype==biotype)
    elif chrom:
        condition = Gene.coord_name==str(chrom)
    elif biotype:
        condition = Gene.biotype==biotype
    else:
        condition = None
    
    query = session.query(Gene)
    if condition is not None:
        query = query.filter(condition)
    
    return query

def get_stable_id_genes_mapping(session):
    """get ensembl genes for a release"""
    genes = session.query(Gene).all()
    ensembl_id_genes = dict([(g.ensembl_id, g) for g in genes])
    return ensembl_id_genes

def get_ranked_expression(session, sample_name, biotype='protein_coding',
        data_path=None, include_target=None, exclude_target=None, rank_by='mean',
        test_run=False):
    """returns all ranked genes from a sample"""
    query = get_gene_expression_query(session, sample_name,
                    data_path=data_path, test_run=test_run)
    
    if biotype:
        query = query.filter(Gene.biotype==biotype)
    
    records = query.all()

    genes = []
    for expressed in records:
        gene = expressed.gene
        gene.Scores = expressed.scores
        genes.append(gene)

    # keep only those genes in the include target gene set if provided
    if include_target is not None:
        include_genes = get_target_genes(session, include_target,
            test_run=test_run)
        if include_genes:
            include_id_set = set([gene.ensembl_id for gene in include_genes])
            final_genes = []
            for gene in genes:  # if it intersects then keep it
                if len(include_id_set.intersection(set([gene.ensembl_id]))) > 0:
                    final_genes.append(gene)
            genes = final_genes

    # keep only those genes not in the exclude target gene set if provided
    if exclude_target is not None:
        exclude_genes = get_target_genes(session, exclude_target,
            test_run=test_run)
        if exclude_genes:
            exclude_id_set = set([gene.ensembl_id for gene in exclude_genes])
            final_genes = []
            for gene in genes: # if it doesn't intersect then keep it
                if len(exclude_id_set.intersection(set([gene.ensembl_id]))) == 0:
                    final_genes.append(gene)
            genes = final_genes

    # set rank
    if rank_by.lower() == 'mean':
        scored = [(g.MeanScore, g) for g in genes]
    elif rank_by.lower() == 'max':
        scored = [(g.MaxScore, g) for g in genes]
    else:
        raise NotImplementedError
    
    # Make sure we get highest first
    scored = reversed(sorted(scored))
    genes = []
    for rank, (score, gene) in enumerate(scored):
        gene.Rank = rank + 1
        genes.append(gene)

    return genes

def get_ranked_expression_diff(session, sample_name, multitest_signif_val=None,
        biotype='protein_coding', data_path=None, include_target=None,
        exclude_target=None, rank_by='mean', test_run=False):
    """returns all ranked genes from a sample difference experiment"""
    query = get_gene_expression_diff_query(session, sample_name,
            data_path=data_path,multitest_signif_val=multitest_signif_val,
            test_run=test_run)

    if biotype:
        query = query.filter(Gene.biotype==biotype)

    records = query.all()

    genes = []
    for expressed_diff in records:
        gene = expressed_diff.gene
        gene.Scores = expressed_diff.fold_changes
        genes.append(gene)

    # keep only those genes in the include target gene set if provided
    if include_target is not None:
        include_genes = get_target_genes(session, include_target,
                test_run=test_run)
        if include_genes:
            include_id_set = set([gene.ensembl_id for gene in include_genes])
            final_genes = []
            for gene in genes:  # if it intersects then keep it
                if len(include_id_set.intersection(set([gene.ensembl_id]))) > 0:
                    final_genes.append(gene)
            genes = final_genes

    # keep only those genes not in the exclude target gene set if provided
    if exclude_target is not None:
        exclude_genes = get_target_genes(session, exclude_target,
            test_run=test_run)
        if exclude_genes:
            exclude_id_set = set([gene.ensembl_id for gene in exclude_genes])
            final_genes = []
            for gene in genes: # if it doesn't intersect then keep it
                if len(exclude_id_set.intersection(set([gene.ensembl_id]))) == 0:
                    final_genes.append(gene)
            genes = final_genes

    # set rank
    if rank_by.lower() == 'mean':
        scored = [(g.MeanScore, g) for g in genes]
    elif rank_by.lower() == 'max':
        scored = [(g.MaxScore, g) for g in genes]
    else:
        raise NotImplementedError

    # Make sure we get highest first
    scored = reversed(sorted(scored))
    genes = []
    for rank, (score, gene) in enumerate(scored):
        gene.Rank = rank + 1
        genes.append(gene)

    return genes

def get_ranked_genes_per_chrom(session, sample_name, chrom, biotype='protein_coding', data_path=None, test_run=False):
    """returns genes from a chromosome"""
    # TODO remove hardcoding for mouse!
    assert chrom in chroms['mouse']
    genes = get_ranked_expression(session, sample_name,
                    biotype=biotype, data_path=data_path, test_run=test_run)
    genes = (g for g in genes if g.coord_name==chrom)
    return tuple(genes)

def get_diff_ranked_genes_per_chrom(session, sample_name, multitest_signif_val, chrom, biotype='protein_coding', data_path=None, test_run=False):
    """returns difference experiment genes from a chromosome"""
    # TODO remove hardcoding for mouse!
    assert chrom in chroms['mouse']
    genes = get_ranked_expression_diff(session, sample_name,
            multitest_signif_val, biotype=biotype, data_path=data_path,
            test_run=test_run)
    genes = (g for g in genes if g.coord_name==chrom)
    return tuple(genes)

def get_target_genes(session, target_gene_sample_name, test_run=False):
    """returns target genes, not ranked"""
    target_sample = _get_sample(session, target_gene_sample_name)
    if not target_sample:
        raise RuntimeError('No target_sample with name %s' % \
                        target_gene_sample_name)
    
    query = session.query(Gene).join(TargetGene).\
            filter(TargetGene.sample_id==target_sample.sample_id)
    return query

def get_expression_diff_genes(session, sample_name, biotype='protein_coding',
        multitest_signif_val=None, test_run=False):
    """returns the expression diff instances"""

    sample = _get_sample(session, sample_name)

    if not sample:
        raise RuntimeError('No sample with name %s' % sample_name)

    query = session.query(ExpressionDiff).\
            filter(ExpressionDiff.sample_id==sample.sample_id)

    query.filter(Gene.biotype==biotype)
    if multitest_signif_val is not None:
        query = query.filter(ExpressionDiff.multitest_signif==multitest_signif_val)

    return query

def get_total_gene_counts(session, sample_name, biotype='protein_coding', data_path=None, test_run=False):
    """docstring for get_total_gene_counts"""
    query = get_gene_expression_query(session, sample_name,
                                        data_path, test_run)
    if biotype:
        query = query.filter(Gene.biotype==biotype)
    
    return query.distinct().count()

def get_total_diff_gene_counts(session, sample_name, biotype='protein_coding', data_path=None, test_run=False):
    """docstring for get_total_gene_counts"""
    query = get_gene_expression_diff_query(session, sample_name,
                                        data_path, test_run)
    if biotype:
        query = query.filter(Gene.biotype==biotype)

    return query.distinct().count()

# TODO: This appears to be an orphaned function
@display_wrap
def diff_expression_study(session, sample_name, data_path, table,
            ensembl_id_label='ENSEMBL', test_run=False,
            run_record=None, ui=None):
    """returns list of genes and why they're not present in an existing
    expression db"""
    if run_record is None:
        run_record = RunRecord()
    # query db for every gene linked to this expression study
    records = get_ranked_expression(session, sample_name,
            data_path=data_path, test_run=test_run)
    # following dead
    # transcript_to_gene = get_transcript_gene_mapping(session)
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
