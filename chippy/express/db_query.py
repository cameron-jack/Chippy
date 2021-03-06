from sqlalchemy.orm.exc import NoResultFound

from chippy.express.db_schema import Chroms, Gene, Exon, TargetGene,\
        Expression, ExpressionDiff, ReferenceFile, Sample, _make_session
from chippy.util.run_record import RunRecord
from chippy.express.util import _one, sample_types

__author__ = 'Gavin Huttley, Cameron Jack'
__copyright__ = 'Copyright 2011-2013, Anuj Pahwa, Gavin Huttley, Cameron Jack'
__credits__ = ['Gavin Huttley', 'Cameron Jack']
__license__ = 'GPL'
__maintainer__ = 'Cameron Jack'
__email__ = 'cameron.jack@anu.edu.au'
__status__ = 'Pre-release'
__version__ = '0.1'

import functools # for use with decorator-style wrappers

### NOTES ###
# This file is structured into a number of function groupings
# 1. Private helper functions
# 2. Public macro functions (complex functionality)
# 3. Public table entry query functions
# 4. Public table entry counting functions
### END NOTES ###

### 1. Private helper functions ###

# Wrapper to ensure _entries functions behave consistently
def _safe_query(query_function):
    """ Wraps a query function to return a list or an empty list """
    @functools.wraps(query_function)
    def safe_query(*args, **kwargs):
        try:
            return query_function(*args, **kwargs)
        except NoResultFound:
            return []
    return safe_query

# Wrapper to ensure _counts functions behave consistently
def _safe_counts(counts_function):
    """ Wraps a count function to return an int or 0 """
    @functools.wraps(counts_function)
    def safe_counts(*args, **kwargs):
        try:
            return counts_function(*args, **kwargs)
        except NameError, AttributeError:
            return 0
    return safe_counts

# Core sample retrieval function
def _get_sample(session, sample_name):
    """ returns a single sample, or None """
    try:
        sample = session.query(Sample).filter(Sample.name==sample_name).one()
    except NoResultFound:
        sample = None
    return sample

### Private get_table_query functions
def _get_reffiles_query(session, reffile_name=None, sample_name=None):
    """ returns a query for reffiles """
    query = session.query(ReferenceFile)
    if reffile_name is not None:
        query = query.filter(ReferenceFile.name==reffile_name)
    if sample_name is not None:
        query = query.join(Sample).filter(Sample.name==sample_name)
    return query

def _get_gene_query(session, biotype='protein_coding', chrom=None, data_path=None):
    """ returns total unique gene entries"""
    query = session.query(Gene)
    if chrom:
        query = query.filter(Gene.chrom==chrom)
    if biotype:
        query = query.filter(Gene.biotype==biotype)
    return query

def _get_exon_query(session, ensembl_id=None, gene_ensembl=None):
    """ returns an exon query object """
    query = session.query(Exon)
    if ensembl_id:
        query = query.filter(Exon.ensembl_id==ensembl_id)
    if gene_ensembl:
        query = query.join(Gene).filter(Gene.ensembl_id==gene_ensembl)
    return query

def _get_expression_query(session, sample_name=None,
        biotype='protein_coding', chrom=None, data_path=None):
    """ Returns expression table query """
    rr = RunRecord('_get_expression_query')
    query = session.query(Expression)
    if sample_name is not None:
        sample = _get_sample(session, sample_name)
        if sample is None:
            rr.dieOnCritical('Unknown sample name', sample_name)
        query = query.filter(Expression.sample_id==sample.sample_id)

    if data_path is not None:
        reffile_id = _one(session.query(ReferenceFile.reffile_id).\
                filter(ReferenceFile.name==data_path))
        if not data_path:
            rr.dieOnCritical('Unknown data path', data_path)
        reffile_id = reffile_id[0]
        query = query.filter(Expression.reffile_id==reffile_id)

        # used to reconstruct the origin of a sample
    query = query.join(Gene)

    if chrom is not None:
        query = query.filter(Gene.chrom==chrom)
    if biotype is not None:
        query = query.filter(Gene.biotype==biotype)
    return query

def _get_diff_query(session, sample_name=None, biotype='protein_coding',
        multitest_signif_val=None, chrom=None, data_path=None):
    """ Returns ExpressionDiff table query """
    rr = RunRecord('_get_diff_query')
    query = session.query(ExpressionDiff)
    if sample_name is not None:
        sample = _get_sample(session, sample_name)
        if not sample:
            rr.dieOnCritical('No sample with name', sample_name)

        query = query.filter(ExpressionDiff.sample_id==sample.sample_id)

    if data_path is not None:
        reffile_id = _one(session.query(ReferenceFile.reffile_id).\
        filter(ReferenceFile.name==data_path))
        if not data_path:
            rr.dieOnCritical('Unknown data path', data_path)
        reffile_id = reffile_id[0]
        query = query.filter(Expression.reffile_id==reffile_id)

    if multitest_signif_val is not None:
        query = query.filter(ExpressionDiff.multitest_signif==\
                                 multitest_signif_val)
    query = query.join(Gene)
    if chrom is not None:
        query = query.filter(Gene.chrom==chrom)
    if biotype:
        query = query.filter(Gene.biotype==biotype)
    return query

def _get_targetgene_query(session, sample_name=None, biotype='protein_coding'):
    """ Returns target_gene records for a given sample """
    rr = RunRecord('get_targets')
    if sample_name is not None:
        sample = _get_sample(session, sample_name)
        if sample is None:
            rr.addError('Using all samples, as no sample matches name',
                    sample_name)
            query = session.query(TargetGene).join(Gene)
        else:
            query = session.query(TargetGene).join(Gene).\
                    filter(TargetGene.sample_id==sample.sample_id)
    else: # get them all
        query = session.query(TargetGene).join(Gene)

    if biotype:
        query = query.filter(Gene.biotype==biotype)
    return query

### 2. Public misc. functions

def make_session(db_path):
    """ Isolates DB type from path """
    session = _make_session('sqlite:///' + db_path)
    return session

def get_all_sample_names(session):
    """returns a list of available choices for samples in the DB"""
    if session is None:
        return ['No connection to DB']
    elif type(session) is str:
        session = make_session(session)
    return [str(s.name) for s in get_sample_entries(session)]

def get_sample_names_descriptions(session):
    """ Returns a dictionary of sample names with descriptions as values """
    if session is None:
        return ['No connection to DB']
    elif type(session) is str:
        session = make_session(session)
    return dict([(str(s.name), str(s.description)) for s in \
            get_sample_entries(session)])

def get_expr_diff_sample_names(session):
    """ returns a list of absolute and differential expression samples """
    if session is None:
        return ['No connection to DB']
    elif type(session) is str:
        session = make_session(session)
    e_names = set([str(e.sample.name) for e in get_expression_entries(session)])
    d_names = set([str(d.sample.name) for d in get_diff_entries(session)])
    return list(e_names.union(d_names))

def get_expr_sample_names(session):
    """ returns a list of available absolute expression samples """
    if session is None:
        return ['No connection to DB']
    elif type(session) is str:
        session = make_session(session)
    expr = get_expression_entries(session)
    return list(set(str(e.sample.name) for e in expr))

def get_diff_sample_names(session):
    """ returns a list of available differential expression samples """
    if session is None:
        return ['No connection to DB']
    elif type(session) is str:
        session = make_session(session)
    diff = get_diff_entries(session)
    return list(set(str(d.sample.name) for d in diff))

def get_target_gene_names(session):
    """ Returns a list of the sample names that are TargetGene entries """
    if session is None:
        return ['No connection to DB']
    elif type(session) is str:
        session = make_session(session)
    targets = get_targetgene_entries(session)
    return list(set(str(t.sample.name) for t in targets))

def get_chroms(session):
    """ return list of chroms from ',' separated string """
    if session is None:
        return ['No connection to DB']
    elif type(session) is str:
        session = make_session(session)
    rr = RunRecord('get_chroms')
    try:
        chroms = session.query(Chroms).one()
        chroms = chroms.chromStr.split(',')
    except NoResultFound:
        chroms = []
        rr.addError('Chroms found', None)
    return chroms

def get_species(session):
    """ returns the species value from Chroms table """
    if session is None:
        return ['No connection to DB']
    elif type(session) is str:
        session = make_session(session)
    try:
        return session.query(Chroms).one().species
    except NoResultFound:
        return None

def get_stable_id_genes_mapping(session):
    """ get ensembl genes as a dict indexed by ensembl_id or empty dict"""
    try:
        genes = session.query(Gene).all()
        ids_genes = dict([(g.ensembl_id, g) for g in genes])
    except NoResultFound:
        ids_genes = {'None':'None'}
    return ids_genes

def get_gene_ids(session, chrom=None, include_targets=None,
        exclude_targets=None):
    """
        Returns a list (or counts) of ensembl_ids restricted by chrom,
        include and exclude target_gene lists.
    """
    # limit by chrom if provided
    genes = get_gene_entries(session, chrom=chrom)
    stable_ids = [g.ensembl_id for g in genes]

    if include_targets is not None:
        include_ids = []
        for target in include_targets:
            genes = get_targetgene_entries(session, sample_name=target)
            for g in genes:
                include_ids.append(g.gene.ensembl_id)
        stable_ids = list(set(stable_ids).\
                intersection(set(include_ids)))

    if exclude_targets is not None:
        exclude_ids = []
        for target in exclude_targets:
            genes = get_targetgene_entries(session, sample_name=target)
            for g in genes:
                exclude_ids.append(g.gene.ensembl_id)
        stable_ids = list(set(stable_ids)-(set(exclude_ids)))

    return stable_ids

def get_genes_by_ranked_expr(session, sample_name, biotype='protein_coding',
        chrom=None, data_path=None, include_targets=None, exclude_targets=None,
        rank_by='mean'):
    """returns all ranked genes from a sample"""
    rr = RunRecord('get_genes_by_ranked_expr')
    records = get_expression_entries(session, sample_name=sample_name,
            biotype=biotype, chrom=chrom, data_path=data_path)

    genes = []
    for expressed in records:
        gene = expressed.gene
        gene.Scores = expressed.scores
        genes.append(gene)

    # keep only those genes in the include target gene set if provided
    if include_targets is not None:
        include_genes = get_targetgene_entries(session, include_targets)
        if len(include_genes) > 0:
            include_gene_ids = set([tg.gene.ensembl_id for tg in include_genes])
            genes = [gene for gene in genes if gene.ensembl_id in \
                    include_gene_ids]

    # keep only those genes not in the exclude target gene set if provided
    if exclude_targets is not None:
        exclude_genes = get_targetgene_entries(session,
                sample_name=exclude_targets)
        if len(exclude_genes) > 0:
            exclude_gene_ids = set([tg.gene.ensembl_id for tg in exclude_genes])
            genes = [gene for gene in genes if gene.ensembl_id not in \
                    exclude_gene_ids]

    # set rank
    if rank_by.lower() == 'mean':
        scored = [(g.MeanScore, g) for g in genes]
    elif rank_by.lower() == 'max':
        scored = [(g.MaxScore, g) for g in genes]
    else:
        rr.dieOnCritical('Ranking method not possible', rank_by.lower())

    # Make sure we get highest first
    scored = reversed(sorted(scored))
    genes = []
    for rank, (score, gene) in enumerate(scored):
        gene.Rank = rank + 1
        genes.append(gene)

    return genes

def get_genes_by_ranked_diff(session, sample_name, multitest_signif_val=None,
        biotype='protein_coding', chrom=None, data_path=None,
        include_targets=None, exclude_targets=None, rank_by='mean'):
    """returns all ranked genes from a sample difference experiment"""
    rr = RunRecord('get_genes_by_ranked_diff')
    records = get_diff_entries(session, sample_name=sample_name,
            biotype=biotype, data_path=data_path, chrom=chrom,
            multitest_signif_val=multitest_signif_val)

    genes = []
    for expressed_diff in records:
        gene = expressed_diff.gene
        gene.Scores = expressed_diff.fold_changes
        genes.append(gene)

    # keep only those genes in the include target gene sets if provided
    if include_targets is not None:
        include_genes = get_targetgene_entries(session, include_targets)
        if len(include_genes) > 0:
            include_gene_ids = set([tg.gene.ensembl_id for tg in include_genes])
            genes = [gene for gene in genes if gene.ensembl_id in\
                    include_gene_ids]

    # keep only those genes not in the exclude target gene sets if provided
    if exclude_targets is not None:
        exclude_genes = get_targetgene_entries(session, exclude_targets)
        if len(exclude_genes) > 0:
            exclude_gene_ids = set([tg.gene.ensembl_id for tg in exclude_genes])
            genes = [gene for gene in genes if gene.ensembl_id not in\
                    exclude_gene_ids]

    # set rank
    if rank_by.lower() == 'mean':
        scored = [(g.MeanScore, g) for g in genes]
    elif rank_by.lower() == 'max':
        scored = [(g.MaxScore, g) for g in genes]
    else:
        rr.dieOnCritical('Ranking method not possible', rank_by.lower())

    # Make sure we get highest first
    scored = reversed(sorted(scored))
    genes = []
    for rank, (score, gene) in enumerate(scored):
        gene.Rank = rank + 1
        genes.append(gene)

    return genes

def drop_sample_records(session, sample_name, test=False):
    """
        Drop the expression/target/expr_diff records, sample name
        and reference file attached to a given sample name.
        Returns True if sample deleted.
    """
    rr = RunRecord('drop_sample_records')
    sample = _get_sample(session, sample_name)
    if sample is not None:
        rr.addInfo('Deleting records for ', sample_name)
        num_removed = session.query(ReferenceFile).filter(ReferenceFile.sample_id==sample.sample_id).delete()
        rr.addInfo('Number of reference file records deleted', num_removed)
        num_removed = session.query(Expression).filter(Expression.sample_id==sample.sample_id).delete()
        rr.addInfo('Number of expression records deleted', num_removed)
        num_removed = session.query(ExpressionDiff).filter(ExpressionDiff.sample_id==sample.sample_id).delete()
        rr.addInfo('Number of expressionDiff records deleted', num_removed)
        num_removed = session.query(TargetGene).filter(TargetGene.sample_id==sample.sample_id).delete()
        rr.addInfo('Number of target gene records deleted', num_removed)
        num_removed = session.query(Sample).filter(Sample.sample_id==sample.sample_id).delete()
        assert num_removed == 1, 'Incorrect number of samples removed' + str(num_removed)
        rr.addInfo('Number of sample records deleted', num_removed)
        try:
            if not test:
                session.commit()
                return True
            else:
                rr.addInfo('session deletes not committed', sample_name)
        except NoResultFound:
            rr.addError('Deletes could not take place', 'Commit failed')
    else:
        rr.addError('No sample found by name', sample_name)
    session.rollback()
    return False


### Public DB table query functions ###
# For querying tables: Sample, ReferenceFile, Gene, Exon, Expression,
# ExpressionDiff, TargetGene (Chrom query functions are in the Public
# misc. functions section)

@_safe_query
def get_sample_entries(session):
    """ Returns all samples, [] on failure """
    return session.query(Sample).all()

@_safe_query
def get_gene_entries(session, biotype='protein_coding', chrom=None, data_path=None):
    """ returns a list of gene entries"""
    query = _get_gene_query(session, biotype=biotype, chrom=chrom,
            data_path=data_path)
    return query.distinct().all()

@_safe_query
def get_expression_entries(session, sample_name=None,
        biotype='protein_coding', chrom=None, data_path=None):
    """ returns a list of Expression entries"""
    query = _get_expression_query(session, sample_name=sample_name,
            biotype=biotype, chrom=chrom, data_path=data_path)
    return query.distinct().all()

@_safe_query
def get_diff_entries(session, sample_name=None,
        biotype='protein_coding', chrom=None,
        multitest_signif_val=None, data_path=None):
    """ Returns a list of unique ExpressionDiff entries """
    query = _get_diff_query(session, sample_name=sample_name,
            biotype=biotype, multitest_signif_val=multitest_signif_val,
            chrom=chrom, data_path=data_path)
    return query.distinct().all()

@_safe_query
def get_targetgene_entries(session, sample_name=None, biotype='protein_coding'):
    """ Returns a list of target_gene records for a given sample """
    query = _get_targetgene_query(session, sample_name, biotype=biotype)
    return query.all()

@_safe_query
def get_reffile_entries(session, reffile_name=None, sample_name=None):
    """ Returns a list of reffiles """
    query = _get_reffiles_query(session, reffile_name=reffile_name,
            sample_name=sample_name)
    return query.all()

@_safe_query
def get_exon_entries(session, ensembl_id=None ,gene_ensembl=None):
    """ returns all exons, or just those for a specific gene stable_id """
    query = _get_exon_query(session, ensembl_id=ensembl_id,
            gene_ensembl=gene_ensembl)
    return query.all()

### Public DB count functions ###
# For counting tables: Sample, ReferenceFile, Gene, Exon, Expression,
# ExpressionDiff, TargetGene (Chrom query functions are in the Public
# misc. functions section)

@_safe_counts
def get_sample_counts(session):
    """ Returns the integer count of samples, 0 on failure """
    return session.query(Sample).count()

@_safe_counts
def get_gene_counts(session, biotype='protein_coding', chrom=None, data_path=None):
    """ returns the number of gene entries"""
    query = _get_gene_query(session, biotype=biotype, chrom=chrom,
            data_path=data_path)
    return query.distinct().count()

@_safe_counts
def get_expression_counts(session, sample_name=None,
        biotype='protein_coding', data_path=None):
    """ returns the number of Expression entries"""
    query = _get_expression_query(session, sample_name=sample_name,
            biotype=biotype, data_path=data_path)
    return query.distinct().count()

@_safe_counts
def get_diff_counts(session, sample_name=None, biotype='protein_coding',
        chrom=None, multi_signif_val=None, data_path=None):
    """ Returns number of unique ExpressionDiff entries """
    query = _get_diff_query(session, sample_name=sample_name,
            biotype=biotype, chrom=chrom,
            multitest_signif_val=multi_signif_val, data_path=data_path)
    return query.distinct().count()

@_safe_counts
def get_targetgene_counts(session, sample_name=None,
        biotype='protein_coding', data_path=None):
    """ Returns target_gene records for a given sample """
    query = _get_targetgene_query(session, sample_name, biotype=biotype)
    return query.count()

@_safe_counts
def get_reffile_counts(session, reffile_name=None, sample_name=None):
    query = _get_reffiles_query(session, reffile_name=reffile_name,
            sample_name=sample_name)
    return query.count()

@_safe_counts
def get_exon_counts(session, ensembl_id=None, gene_ensembl=None):
    """ returns counts of all exons, or just those for a specific gene
            stable_id """
    query = _get_exon_query(session, ensembl_id=ensembl_id,
            gene_ensembl=gene_ensembl)
    return query.count()

### END of db_query

