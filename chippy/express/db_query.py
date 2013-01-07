from sqlalchemy import and_
from sqlalchemy.orm import contains_eager, joinedload
from sqlalchemy.orm.exc import NoResultFound

from chippy.express.db_schema import Chroms, Gene, Exon, TargetGene,\
        Expression, ExpressionDiff, ReferenceFile, Sample, _make_session
from chippy.util.run_record import RunRecord
from chippy.express.util import single_gene, _one

__author__ = "Gavin Huttley, Cameron Jack"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley, Cameron Jack"
__credits__ = ["Gavin Huttley, Cameron Jack"]
__license__ = "GPL"
__maintainer__ = "Cameron Jack"
__email__ = "cameron.jack@anu.edu.au"
__status__ = "Pre-release"
__version__ = '0.1'

### Public interface ###

def make_session(db_path):
    session = _make_session(db_path)
    return session

def get_sample_choices(session):
    """returns the available choices for samples in the DB"""
    samples = [str(s.name) + ' : ' + str(s.description)\
               for s in get_samples(session)]
    # These are valid null samples
    samples.append('None')
    session.close()
    return samples

def get_chroms(session, species):
    try:
        chroms = session.query(Chroms).filter(Chroms.species==species).one()
        chroms = chroms.chromStr.split(',')
    except NoResultFound:
        print 'No chroms found for species:', species
        chroms = []
    return chroms

def get_samples(session):
    """returns all samples"""
    try:
        samples = session.query(Sample).all()
    except NoResultFound:
        samples = []
    return samples

def get_samples_count(session):
    """ returns the integer count for total number of samples """
    try:
        count = get_samples(session).count()
    except NoResultFound:
        count = 0
    return count

def get_target_samples(session):
    """returns samples for target genes"""
    try:
        query = session.query(Sample).join(TargetGene).all()
    except NoResultFound:
        query = []
    return query

def get_genes(session, chrom=None, biotype='protein_coding',
        stable_ids=None):
    """returns the Gene's for the indicated release, chrom,
    biotype OR ensembl stable ID.
    
    Note: if stable_ids provided, all arguments other are ignored.
    """
    if type(stable_ids) == str:
        stable_ids = [stable_ids]
    
    if stable_ids:
        condition = Gene.ensembl_id.in_(stable_ids)
    elif chrom and biotype:
        condition = and_(Gene.coord_name==str(chrom),
                Gene.biotype==biotype)
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
    """get ensembl genes as a dict indexed by ensembl_id"""
    try:
        genes = session.query(Gene).all()
        ids_genes = dict([(g.ensembl_id, g) for g in genes])
    except NoResultFound:
        ids_genes = {'None':'None'}
    return ids_genes

def get_ranked_abs_expr_genes(session, sample_name,
        biotype='protein_coding', data_path=None, include_target=None,
        exclude_target=None, rank_by='mean', test_run=False):
    """returns all ranked genes from a sample"""
    query = _get_gene_expression_query(session, sample_name,
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
        try:
            include_gene_ids = set([gene.ensembl_id for gene in \
                    get_target_genes(session, include_target, \
                    test_run=test_run)])
        except NoResultFound:
            include_gene_ids = None
        if include_gene_ids:
            genes = [gene for gene in genes if gene.ensembl_id in \
                    include_gene_ids]

    # keep only those genes not in the exclude target gene set if provided
    if exclude_target is not None:
        try:
            exclude_gene_ids = set([gene.ensembl_id for gene in \
                    get_target_genes(session, exclude_target,
                    test_run=test_run)])
        except NoResultFound:
            exclude_gene_ids = None
        if exclude_gene_ids:
            genes = [gene for gene in genes if gene.ensembl_id in \
                    exclude_gene_ids]

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

def get_ranked_diff_expr_genes(session, sample_name, multitest_signif_val=None,
        biotype='protein_coding', data_path=None, include_target=None,
        exclude_target=None, rank_by='mean', test_run=False):
    """returns all ranked genes from a sample difference experiment"""
    query = _get_gene_expression_diff_query(session, sample_name,
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

def get_ranked_genes_per_chrom(session, sample_name, species, chrom,
        biotype='protein_coding', data_path=None, test_run=False):
    """returns genes from a chromosome"""
    assert chrom in get_chroms(session, species)
    genes = get_ranked_abs_expr_genes(session, sample_name,
            biotype=biotype, data_path=data_path, test_run=test_run)
    genes = (g for g in genes if g.coord_name==chrom)
    return tuple(genes)

def get_diff_ranked_genes_per_chrom(session, sample_name,
        multitest_signif_val, species, chrom, biotype='protein_coding',
        data_path=None, test_run=False):
    """returns difference experiment genes from a chromosome"""
    assert chrom in get_chroms(session, species)
    genes = get_ranked_diff_expr_genes(session, sample_name,
            multitest_signif_val, biotype=biotype, data_path=data_path,
            test_run=test_run)
    genes = (g for g in genes if g.coord_name==chrom)
    return tuple(genes)

def get_target_genes(session, target_gene_sample_name, test_run=False):
    """returns target genes, not ranked"""
    target_sample = _get_sample(session, target_gene_sample_name)
    if not target_sample:
        raise RuntimeError('No target_sample with name ' +\
                target_gene_sample_name)
    query = session.query(Gene).join(TargetGene).\
            filter(TargetGene.sample_id==target_sample.sample_id)
    return query

def get_expr_entries(session, sample_name=None, biotype='protein_coding',
        count_only=False, test_run=False):
    """ Returns expression records for a given sample """
    if sample_name:
        sample = _get_sample(session, sample_name)
        if not sample:
            raise RuntimeError('No sample with name ' + sample_name)

        query = session.query(Expression.expression_id).\
                filter(Expression.sample_id==sample.sample_id)
        #if biotype:
        #    query = query.filter(Gene.biotype==biotype)

        if count_only:
            try:
                count = query.count()
            except NoResultFound:
                count = 0
            return count
        else:
            return query.all()
    else:
        query = session.query(Expression.expression_id)
        #if biotype:
        #    query = query.filter(Gene.biotype==biotype)

        if count_only:
            try:
                count = query.count()
            except NoResultFound:
                count = 0
            return count
        else:
            return query.all()

def get_diff_entries(session, sample_name=None, biotype='protein_coding',
        multitest_signif_val=None, count_only=False, test_run=False):
    """ Returns expression_diff records for a given sample """

    if sample_name:
        sample = _get_sample(session, sample_name)
        if not sample:
            raise RuntimeError('No sample with name ' + sample_name)

        query = session.query(ExpressionDiff.expression_diff_id).\
                filter(ExpressionDiff.sample_id==sample.sample_id)
        #if biotype:
        #    query = query.filter(Gene.biotype==biotype)
        if multitest_signif_val is not None:
            query = query.filter(ExpressionDiff.multitest_signif==\
                    multitest_signif_val)
        if count_only:
            try:
                count = query.count()
            except NoResultFound:
                count = 0
            return count
        else:
            return query.all()
    else:
        query = session.query(ExpressionDiff.expression_diff_id)
        #if biotype:
        #    query = query.filter(Gene.biotype==biotype)
        if multitest_signif_val is not None:
            query = query.filter(ExpressionDiff.multitest_signif==\
                    multitest_signif_val)
        if count_only:
            try:
                count = query.count()
            except NoResultFound:
                count = 0
            return count
        else:
            return query.all()

def get_target_entries(session, sample_name=None, count_only=False, test_run=False):
    """ Returns target_gene records for a given sample """
    if sample_name:
        sample = _get_sample(session, sample_name)
        if not sample:
            raise RuntimeError('No sample with name ' + sample_name)
        query = session.query(TargetGene.target_gene_id).\
                filter(TargetGene.sample_id==sample.sample_id)
        if count_only:
            try:
                count = query.count()
            except NoResultFound:
                count = 0
            return count
        else:
            return query.all()

    else: # get them all
        query = session.query(TargetGene.target_gene_id)
    if count_only:
        try:
            count = query.count()
        except NoResultFound:
            count = 0
        return count
    else:
        return query.all()

### Public total gene number count DB functions ###

def get_total_gene_count(session, sample_name,
        biotype='protein_coding', data_path=None, test_run=False):
    """ returns the number of gene entries"""
    try:
        query = _get_gene_expression_query(session, sample_name,
                data_path, test_run)
        if biotype:
            query = query.filter(Gene.biotype==biotype)
        count = query.distinct().count()
    except NoResultFound:
        count = 0
    return count

def get_total_diff_gene_count(session, sample_name,
        biotype='protein_coding', data_path=None, test_run=False):
    """ returns the number of diff gene entries """
    try:
        query = _get_gene_expression_diff_query(session, sample_name,
                data_path, test_run)
        if biotype:
            query = query.filter(Gene.biotype==biotype)
        count = query.distinct().count()
    except NoResultFound:
        count = 0
    return count

def get_total_target_gene_count(session, sample_name,
        biotype='protein_coding', data_path=None, test_run=False):
    """ returns the number of target gene entries """
    try:
        query = session.query(TargetGene)
        if biotype:
            query = query.filter(Gene.biotype==biotype)
        count = query.distinct().count()
    except NoResultFound:
        count = 0
    return count

### Private helper functions ###

def _get_sample(session, sample_name):
    """ returns a single sample, or None """
    try:
        sample = session.query(Sample).filter(Sample.name==sample_name).one()
    except NoResultFound:
        sample = None
    return sample

def _get_gene_expression_query(session, sample_name, data_path=None,
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

def _get_gene_expression_diff_query(session, sample_name, data_path=None,
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