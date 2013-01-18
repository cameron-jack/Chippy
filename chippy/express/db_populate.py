import warnings
warnings.filterwarnings('ignore', 'Not using MPI as mpi4py not found')

import datetime, sys
sys.path.extend(['..'])
from sqlalchemy import create_engine, and_
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm.exc import NoResultFound

from cogent.db.ensembl import HostAccount, Genome
from cogent.util.progress_display import display_wrap

from chippy.express.db_schema import Chroms, Gene, Exon, \
            TargetGene, Expression, ExpressionDiff, ReferenceFile, Sample
from chippy.express.db_query import  get_stable_id_genes_mapping, \
        get_expr_entries, get_diff_entries, get_target_entries
from chippy.express.util import sample_types, _one
from chippy.util.run_record import RunRecord
from chippy.util.definition import LOG_DEBUG, LOG_INFO, LOG_WARNING, \
    LOG_ERROR, LOG_CRITICAL


__author__ = "Gavin Huttley, Cameron Jack"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley, Cameron Jack"
__credits__ = ["Gavin Huttley, Cameron Jack"]
__license__ = "GPL"
__maintainer__ = "Cameron Jack"
__email__ = "cameron.jack@anu.edu.au"
__status__ = "Pre-release"
__version__ = '0.1'

now = datetime.datetime.now()
today = datetime.date(now.year, now.month, now.day)

def successful_commit(session, data, debug=False):
    """returns True if successfully added data"""
    if not data:
        return False
    
    try:
        data[0]
    except TypeError:
        data = [data]
    
    if debug:
        print data

    try:
        session.add_all(data)
        session.commit()
    except IntegrityError, msg:
        session.rollback()
        if debug:
            print msg
        return False
    return True

def add_chroms(session, species, chromlist):
    """ adds a species and its list of chroms to the DB.
    This interface is for external use """
    chromsList = [str(chrom) for chrom in chromlist]
    chroms = Chroms(species, chromsList)
    success = successful_commit(session, chroms)
    if not success:
        return False
    return True

def add_ensembl_gene_data(session, species, ensembl_release, account=None,
        rr=RunRecord(), debug=False):
    """add Ensembl genes and their transcripts to the db session"""
    
    genome = Genome(species, Release=ensembl_release, account=account)
    skip = set(['processed_transcript', 'pseudogene'])
    biotypes = [b for b in genome.getDistinct('BioType') if b not in skip]
    
    data = []
    unique_gene_ids = set()
    unique_exon_ids = set()
    chromSet = set()
    n = 0
    total_objects = 0
    for biotype in biotypes:
        for gene in genome.getGenesMatching(BioType=biotype):
            # gene.Location.CoordName is the chromosome name
            chromSet.add(gene.Location.CoordName)

            if gene.StableId not in unique_gene_ids:

                db_gene = Gene(ensembl_id=gene.StableId,
                        symbol=gene.Symbol, biotype=gene.BioType,
                        description=gene.Description, status=gene.Status,
                        chrom=gene.Location.CoordName,
                        start=gene.Location.Start, end=gene.Location.End,
                        strand=gene.Location.Strand)
                total_objects += 1
                unique_gene_ids.add(gene.StableId)
                data.append(db_gene)
            else:
                rr.addWarning('add_ensembl_gene_data', 'Duplicate gene',
                        gene.StableId)

            for exon in gene.CanonicalTranscript.Exons:
                if exon.StableId not in unique_exon_ids:
                    db_exon = Exon(exon.StableId, exon.Rank,
                            exon.Location.Start, exon.Location.End)
                    db_exon.gene = db_gene
                    unique_exon_ids.add(exon.StableId)
                    data.append(db_exon)
                    total_objects += 1
                else:
                    rr.addWarning('add_ensembl_gene_data', 'Duplicate exon',
                            exon.StableId)
            n += 1
            if n % 100 == 0:
                print 'Genes processed=%s; Db objects created=%d' % (n,
                        total_objects)
                if debug:
                    session.add_all(data)
                    session.commit()
                    return
    rr.addInfo('add_ensembl_gene_data', 'Instantiating chromosomes',
            chromSet)
    chroms = Chroms(species, chromSet)
    data.append(chroms)

    rr.addInfo('add_ensembl_gene_data', 'Writing objects into db',
            total_objects)
    session.add_all(data)
    session.commit()
    return chroms, rr

def add_sample(session, name, description, rr=RunRecord()):
    """add a basic sample"""
    sample = Sample(name, description)
    if not successful_commit(session, sample):
        rr.addMessage('add_samples', LOG_INFO,
                'Sample already exists in db', name)
        return False, rr
    else:
        rr.addMessage('add_samples', LOG_INFO,
                'Sample created in db', name)
        return True, rr

@display_wrap
def add_expression_study(session, sample_name, data_path, table,
        probeset_label='probeset', ensembl_id_label='ENSEMBL',
        expression_label='exp', rr=RunRecord(),
        ui=None):
    """adds Expression instances into the database from table
    
    Arguments:
        - sample_name: the sample name to connect expression instances to
        - data_path: the reference file path
        - table: the actual expression data table
        - probeset_label: label of the column containing probset id
        - ensembl_id_label: label of the column containing Ensembl Stable IDs
        - expression_label: label of the column containing absolute measure
          of expression
    """

    data = []
    sample = _one(session.query(Sample).filter_by(name=sample_name))
    if not sample:
        session.rollback()
        raise RuntimeError('error querying for a sample')
    
    reffile = session.query(ReferenceFile).filter_by(name=data_path).all()
    if len(reffile) == 0:
        reffile = ReferenceFile(data_path, today)
        reffile.sample = sample
        data.append(reffile)
    else:
        reffile = reffile[0]
    
    # check we haven't already added expression data from this file, for
    # this sample
    records = session.query(Expression).filter(
            and_(Expression.reffile_id==reffile.reffile_id,
                Expression.sample_id==sample.sample_id)).all()
    
    if len(records) > 0:
        rr.addMessage('add_expression_study',
            LOG_WARNING,
            'Already added this data for this sample / file combo',
            (sample.name, data_path))
        return False, rr
    
    # get all gene ID data for the specified Ensembl release
    ensembl_genes = get_stable_id_genes_mapping(session)
    if not successful_commit(session, data):
        pass
    
    data = []
    unknown_ids = 0
    total = 0
    for record in ui.series(table, noun='Adding expression instances'):
        ensembl_id = record[ensembl_id_label]
        try:
            gene = ensembl_genes[ensembl_id]
        except KeyError:
            unknown_ids += 1
            continue
        
        probeset = record[probeset_label]
        scores = record[expression_label]
        expressed = Expression(probeset, scores)
        expressed.reffile_id = reffile.reffile_id
        expressed.sample = sample
        expressed.gene = gene
        data.append(expressed)
        total += 1

    try:
        session.add_all(data)
        session.commit()
    except IntegrityError:
        session.rollback()
        return False, rr

    if unknown_ids:
        rr.addMessage('add_expression_study',
                LOG_ERROR, 'Number of unknown gene Ensembl IDs', unknown_ids)
    rr.addMessage('add_expression_study',
        LOG_INFO, 'Total genes', total)
    return True, rr

@display_wrap
def add_expression_diff_study(session, sample_name, data_path, table,
            ref_a_path, ref_b_path,
            probeset_label='probeset',
            ensembl_id_label='ENSEMBL', expression_label='exp',
            prob_label='rawp', sig_label='sig',
            rr=RunRecord(), ui=None):
    """adds Expression instances into the database from table
    
    Arguments:
        - sample_name: a description of this sample, e.g. "M-S"
        - data_path: the reference file path
        - table: the actual expression data table
        - ref_a_path, ref_b_path: the difference file contains measurements
          between file path a and file path b. This order is CRITICAL as it
          affects the inferences concerning gene expression being up / down
          in a sample.
        - probeset_label: label of the column containing probset id
        - ensembl_id_label: label of the column containing Ensembl Stable IDs
        - expression_label: label of the column containing absolute measure
          of expression
        - prob_label: label of the column containing raw probability of
          difference in gene expression between samples A/B
        - sig_label: label of the column classifying probabilities as
          significant after multiple test correction. 1 means up in A relative
          to B, -1 means down in A relative to B, 0 means no difference.
        - rr: a RunRecord instance for tracking messages. If not
          provided, one is created and returned.
    """
    
    sample = _one(session.query(Sample).filter_by(name=sample_name))
    if not sample:
        session.rollback()
        raise RuntimeError('error querying for a sample')
    
    ref_a = _one(session.query(ReferenceFile).filter_by(name=ref_a_path))
    if not ref_a:
        rr.addMessage('add_expression_diff_study', LOG_WARNING,
                'Could not find a record for ref_a %s' % ref_a_path,
                str(session.query(ReferenceFile).\
                filter_by(name=ref_a_path).all()))
    
    ref_b = _one(session.query(ReferenceFile).filter_by(name=ref_b_path))
    if not ref_b:
        rr.addMessage('add_expression_diff_study', LOG_WARNING,
                'Could not find a record for ref_b %s' % ref_b_path,
                str(session.query(ReferenceFile).\
                filter_by(name=ref_b_path).all()))
    
    if not ref_a or not ref_b:
        rr.display()
        raise RuntimeError('Reference files not added yet?')
    
    data = []
    reffile = session.query(ReferenceFile).filter_by(name=data_path).all()
    if len(reffile) == 0:
        reffile = ReferenceFile(data_path, today, ref_a_name=ref_a_path,
                ref_b_name=ref_b_path)
        reffile.sample = sample
        data.append(reffile)
    else:
        reffile = reffile[0]

    # check we haven't already added expression data from this file, for
    # this sample
    records = session.query(ExpressionDiff).filter(
            and_(ExpressionDiff.reffile_id==reffile.reffile_id,
            ExpressionDiff.sample_id==sample.sample_id)).all()

    if len(records) > 0:
        rr.addMessage('add_expression_diff_study', LOG_WARNING,
                'Already added this data for this sample / file combo',
                (sample.name, data_path))
        return rr
    
    if not successful_commit(session, data):
        session.rollback()
    
    ensembl_genes = get_stable_id_genes_mapping(session)
    data = []
    unknown_ids = 0
    total = 0
    signif_up_total = 0
    signif_down_total = 0
    for record in ui.series(table, noun='Adding expression diffs'):
        ensembl_id = record[ensembl_id_label]
        try:
            gene = ensembl_genes[ensembl_id]
        except KeyError:
            unknown_ids += 1
            continue
        
        probeset = record[probeset_label]
        fold_change = record[expression_label]
        prob = float(record[prob_label])
        signif = int(record[sig_label])
        diff = ExpressionDiff(probeset, fold_changes=fold_change,
                prob=prob, signif=signif)
        diff.reffile_id = reffile.reffile_id
        diff.sample = sample
        diff.gene = gene
        data.append(diff)
        total += 1
        if signif == 1:
            signif_up_total += 1
        elif signif == -1:
            signif_down_total += 1
    
    session.add_all(data)
    session.commit()

    if unknown_ids:
        rr.addMessage('add_expression_diff_study',
                LOG_ERROR, 'Number of unknown gene Ensembl IDs',
                unknown_ids)
    rr.addMessage('add_expression_diff_study',
            LOG_INFO, 'Total significantly up-regulated genes',
            signif_up_total)
    rr.addMessage('add_expression_diff_study',
            LOG_INFO, 'Total significantly down-regulated genes',
            signif_down_total)
    rr.addMessage('add_expression_diff_study',
            LOG_INFO, 'Total genes', total)
    
    return rr

def _chunk_id_list(id_list, n):
    for i in xrange(0, len(id_list), n):
        yield id_list[i:i+n]

def add_target_genes(session, sample_name, data_path, table,
        ensembl_id_label='ENSEMBL', rr=RunRecord()):
    """adds Expression instances into the database from table

    Arguments:
        - data_path: the reference file path
        - table: the actual expression data table
        - ensembl_id_label: label of the column containing Ensembl Stable IDs
    """
    data = []
    sample = _one(session.query(Sample).filter_by(name=sample_name))
    if not sample:
        session.rollback()
        raise RuntimeError('error querying for a sample')
    
    reffile = session.query(ReferenceFile).filter_by(name=data_path).all()
    if len(reffile) == 0:
        reffile = ReferenceFile(data_path, today)
        reffile.sample = sample
        data.append(reffile)
    else: # Don't overwrite anything, exit instead
        rr.addWarning('add_target_genes', 'File already loaded', data_path)
        rr.display()
        sys.exit(-1)
        #reffile = reffile[0]
    
    ensembl_ids = table.getRawData(ensembl_id_label)

    for id_chunk in _chunk_id_list(ensembl_ids, 100):
        genes = session.query(Gene).filter(Gene.ensembl_id.in_(id_chunk)).all()
        for gene in genes:
            target = TargetGene()
            target.gene = gene
            target.reference_file = reffile
            target.sample = sample
            data.append(target)

    rr.addMessage('add_target_genes', LOG_INFO,
            'Added target genes from', data_path)
    rr.addMessage('add_target_genes', LOG_INFO,
            'No. genes added', len(data))
    session.add_all(data)
    session.commit()
    return rr

def check_existing_data(session, sample_name):
    """ Check if data is already loaded for a given sample name.
    Return the existing sample type and data set size (or None/None) """
    existing_data = get_expr_entries(session, sample_name, count_only=True)
    if existing_data > 0:
        existing_type = sample_types['exp_absolute']
        return existing_data, existing_type

    existing_data = get_diff_entries(session, sample_name, count_only=True)
    if existing_data > 0:
            existing_type = sample_types['exp_absolute']
            return existing_data, existing_type

    existing_data = get_target_entries(session, sample_name, count_only=True)
    if existing_data > 0:
        existing_type = sample_types['target_genes']
        return existing_data, existing_type

    return 0, None

def add_data(session, name, description, path, expr_table,
        gene_id_heading='gene', probeset_heading='probeset',
        expr_heading='exp', sample_type=sample_types['exp_absolute'],
        reffile1=None, reffile2=None, rr=RunRecord()):
    """ A unified interface for adding data to the DB """

    success, rr = add_sample(session, name, description, rr=rr)
    if not success:
        # Check if any sample exists without data
        existing_data, existing_type = check_existing_data(session, name)
        if existing_data > 0:
            rr.addError('add_data', name + ' already has data loaded',
                    existing_data)
            rr.addError('add_data', 'data of type', existing_type)
            return False, rr
        else:
            rr.addInfo('add_data', 'now loading data for existing '\
                    'sample', name)

    # either sample was created or existed with no data, so load data now
    if sample_types[sample_type] == sample_types['exp_absolute']:
        success, rr = add_expression_study(session, name, path, expr_table,
                probeset_label=probeset_heading,
                ensembl_id_label=gene_id_heading,
                expression_label=expr_heading,
                rr=rr)
    elif sample_types[sample_type] == sample_types['exp_diff']:
        # diff between two files, check we got the related files
        assert reffile1 is not None and reffile2 is not None,\
        'To enter differences in gene expression you must specify the 2'\
        'files that contain the absolute measures.'
        rr = add_expression_diff_study(session, name, path,
                expr_table, reffile1, reffile2,
                probeset_label=probeset_heading,
                ensembl_id_label=gene_id_heading,
                expression_label=expr_heading,
                prob_label='rawp', sig_label='sig', rr=rr)
    elif sample_types[sample_type] == sample_types['target_genes']:
        rr = add_target_genes(session, name, path, expr_table,
                ensembl_id_label=gene_id_heading, rr=rr)
    else:
        rr.display()
        raise RuntimeError('Unknown sample type')

    return success, rr

def create_dummy_flat_expr(session, rr=RunRecord()):
    """ create flat and spread dummy data """
    genes_dict = get_stable_id_genes_mapping(session)

    # flat expression dummy
    expr_table_rows = []
    for i, gene_id in enumerate(genes_dict):
        #expr_table_rows.append([gene_id, 'P'+str(i), 1])
        expr_table_row = {}
        expr_table_row['gene'] = gene_id
        expr_table_row['probeset'] = 'P'+str(i)
        expr_table_row['exp'] = 1
        expr_table_rows.append(expr_table_row)

    success, rr = add_data(session, 'dummy_flat',
            'each gene has expression score of 1',
            'dummy_flat_expr.fake', expr_table_rows, gene_id_heading='gene',
            probeset_heading='probeset', expr_heading='exp',
            sample_type='exp_absolute',
            reffile1=None, reffile2=None, rr=rr)

    return success, rr

def create_dummy_spread_expr(session, rr=RunRecord()):
    """ create flat and spread dummy data """
    genes_dict = get_stable_id_genes_mapping(session)

    # spread expression dummy
    expr_table_rows = []
    for i, gene_id in enumerate(genes_dict):
        expr_table_row = {}
        expr_table_row['gene'] = gene_id
        expr_table_row['probeset'] = 'P'+str(i)
        expr_table_row['exp'] = i
        expr_table_rows.append(expr_table_row)

    success, rr = add_data(session, 'dummy_spread',
            'each gene has unique expression score',
            'dummy_spread_expr.fake', expr_table_rows,
            gene_id_heading='gene',
            probeset_heading='probeset', expr_heading='exp',
            sample_type='exp_absolute',
            reffile1=None, reffile2=None, rr=rr)

    return success, rr

