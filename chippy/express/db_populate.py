import datetime
from sqlalchemy import create_engine, and_
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm.exc import NoResultFound

from cogent import LoadTable
from cogent.db.ensembl import HostAccount, Genome
from cogent.util.progress_display import display_wrap

from chippy.express.db_schema import Gene, Exon, \
            ExternalGene, Expression, ExpressionDiff, ReferenceFile, Sample, \
            Session, Base, make_session
from chippy.express.db_query import get_ensembl_id_transcript_mapping, \
        get_transcript_gene_mapping
from chippy.ref.util import chroms
from chippy.express.util import single_gene, _one
from chippy.parse.r_dump import RDumpToTable
from chippy.util.run_record import RunRecord
from chippy.util.definition import LOG_DEBUG, LOG_INFO, LOG_WARNING, \
    LOG_ERROR, LOG_CRITICAL


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"
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
    
    session.add_all(data)
    try:
        session.commit()
    except IntegrityError, msg:
        if debug:
            print msg
        return False
    return True

def add_ensembl_gene_data(session, species, ensembl_release, account=None, debug=False):
    """add Ensembl genes and their transcripts to the db session"""
    
    species_chroms = chroms[species]
    
    genome = Genome(species, Release=ensembl_release, account=account)
    skip = set(['processed_transcript', 'pseudogene'])
    biotypes = [b for b in genome.getDistinct('BioType') if b not in skip]
    
    data = []
    n = 0
    total_objects = 0
    for biotype in biotypes:
        for gene in genome.getGenesMatching(BioType=biotype):
            if gene.Location.CoordName not in species_chroms:
                continue
            
            db_gene = Gene(ensembl_id=gene.StableId, symbol=gene.Symbol,
                biotype=gene.BioType,
                description=gene.Description, status=gene.Status,
                coord_name=gene.Location.CoordName,
                start=gene.Location.Start, end=gene.Location.End,
                strand=gene.Location.Strand,
                ensembl_release=ensembl_release)
            total_objects += 1
            data.append(db_gene)
            for exon in gene.CanonicalTranscript.Exons:
                db_exon = Exon(exon.StableId, exon.Rank,
                        exon.Location.Start, exon.Location.End,
                        ensembl_release)
                db_exon.gene = db_gene
                data.append(db_exon)
                total_objects += 1
                n += 1
            
            if n % 100 == 0:
                print 'Genes processed=%s; Db objects created=%d' % (n,
                                                            total_objects)
                if debug:
                    session.add_all(data)
                    session.commit()
                    return
                
            
    
    print 'Writing objects into db'
    session.add_all(data)
    session.commit()

def add_samples(session, names_descriptions, run_record=None):
    """add basic sample type"""
    if run_record is None:
        run_record = RunRecord()
    
    for name, description in names_descriptions:
        sample = Sample(name, description)
        if not successful_commit(session, sample):
            run_record.addMessage('add_samples',
                LOG_INFO, 'Sample already exists in db', name)
            
            session.rollback()
    
    return run_record


@display_wrap
def add_expression_study(session, sample_name, data_path, table,
        ensembl_release='58', probeset_label='probeset',
        ensembl_id_label='ENSEMBL', expression_label='exp', run_record=None,
        ui=None):
    """adds Expression instances into the database from table
    
    Arguments:
        - sample_name: the sample name to connect expression instances to
        - data_path: the reference file path
        - table: the actual expression data table
        - ensembl_release: the Ensembl release
        - probeset_label: label of the column containing probset id
        - ensembl_id_label: label of the column containing Ensembl Stable IDs
        - expression_label: label of the column containing absolute measure
          of expression
    """
    if run_record is None:
        run_record = RunRecord()
    
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
    
    # check we haven't already added expression data from this file
    records = session.query(Expression).filter_by(reference_file=reffile).all()
    if len(records) > 0:
        run_record.addMessage('add_expression_study',
            LOG_WARNING, 'Already added this expression file', data_path)
        return run_record
    
    # get all gene ID data for the specified Ensembl release
    ensembl_to_transcript = get_ensembl_id_transcript_mapping(session,
                                                            ensembl_release)
    transcript_to_gene = get_transcript_gene_mapping(session, ensembl_release)
    if not successful_commit(session, data):
        pass
    
    data = []
    # order the table in descending order of expression
    table = table.sorted(columns=expression_label, reverse=expression_label)
    rank = 0
    failed = []
    all_genes = set()
    probeset_to_many = 0
    many_probeset = 0
    for record in ui.series(table, noun='Adding expression instances'):
        transcript_ids = record[ensembl_id_label]
        # only keep data from probesets that map to a single gene
        gene_id = single_gene(transcript_to_gene, transcript_ids)
        if gene_id is None:
            probeset_to_many += 1
            continue
        
        all_genes.update([gene_id])
        probeset = record[probeset_label]
        expressed = Expression(probeset, record[expression_label], rank)
        expressed.reffile_id = reffile.reffile_id
        expressed.sample = sample
        data = [expressed]
        
        for ensembl_id in transcript_ids:
            transcript = ensembl_to_transcript[ensembl_id]
            a = Association()
            a.expressed = expressed
            a.transcript_id = transcript.transcript_id
            transcript.expressions.append(a)
            data.extend([a, transcript])
        
        if not successful_commit(session, data):
            many_probeset += 1
            session.rollback()
            continue
        
        rank += 1
    
    session.commit()
    
    run_record.addMessage('add_expression_study',
        LOG_INFO, 'Probset mapped to multiple loci', probeset_to_many)
    run_record.addMessage('add_expression_study',
        LOG_ERROR, 'other probeset error', many_probeset)
    run_record.addMessage('add_expression_study',
        LOG_INFO, 'Total probesets added to db', rank)
    run_record.addMessage('add_expression_study',
        LOG_INFO, 'Total genes uniquely referenced', len(all_genes))
    return run_record

@display_wrap
def add_expression_diff_study(session, data_path, table,
            ref_a_path, ref_b_path,
            ensembl_release='58', probeset_label='probeset',
            ensembl_id_label='ENSEMBL', expression_label='exp',
            prob_label='rawp', sig_label='sig',
            run_record=None, ui=None):
    """adds Expression instances into the database from table
    
    Arguments:
        - data_path: the reference file path
        - table: the actual expression data table
        - ref_a_path, ref_b_path: the difference file contains measurements
          between file path a and file path b. This order is CRITICAL as it
          affects the inferences concerning gene expression being up / down
          in a sample.
        - ensembl_release: the Ensembl release
        - probeset_label: label of the column containing probset id
        - ensembl_id_label: label of the column containing Ensembl Stable IDs
        - expression_label: label of the column containing absolute measure
          of expression
        - prob_label: label of the column containing raw probability of
          difference in gene expression between samples A/B
        - sig_label: label of the column classifying probabilities as
          significant after multiple test correction. 1 means up in A relative
          to B, -1 means down in A relative to B, 0 means no difference.
        - run_record: a RunRecord instance for tracking messages. If not
          provided, one is created and returned.
    """
    
    if run_record is None:
        run_record = RunRecord()
    
    ref_a = _one(session.query(ReferenceFile).filter_by(name=ref_a_path))
    if not ref_a:
        run_record.addMessage('add_expression_diff_study',
            LOG_WARNING, 'Could not find a record for ref_a %s' % ref_a_path,
            str(session.query(ReferenceFile).filter_by(name=ref_a_path).all()))
    
    ref_b = _one(session.query(ReferenceFile).filter_by(name=ref_b_path))
    if not ref_b:
        run_record.addMessage('add_expression_diff_study',
            LOG_WARNING, 'Could not find a record for ref_b %s' % ref_b_path,
            str(session.query(ReferenceFile).filter_by(name=ref_b_path).all()))
    
    if not ref_a or not ref_b:
        run_record.display()
        raise RuntimeError('Reference files not added yet?')
    
    sample_a = ref_a.sample
    sample_b = ref_b.sample
    
    data = []
    reffile = session.query(ReferenceFile).filter_by(name=data_path).all()
    if len(reffile) == 0:
        reffile = ReferenceFile(data_path, today)
        data.append(reffile)
    else:
        reffile = reffile[0]
    
    if not successful_commit(session, data):
        session.rollback()
    
    
    # get all transcript to gene ID data for the specified Ensembl release
    ensembl_to_transcript = get_ensembl_id_transcript_mapping(session,
                                                            ensembl_release)
    transcript_to_gene = get_transcript_gene_mapping(session, ensembl_release)
    
    table = table.sorted(columns=expression_label, reverse=expression_label)
    num = 0
    all_genes = set()
    probeset_to_many = 0
    many_probeset = 0
    for record in ui.series(table, noun='Adding expression diffs'):
        transcript_ids = record[ensembl_id_label]
        gene_id = single_gene(transcript_to_gene, transcript_ids)
        if gene_id is None:
            probeset_to_many += 1
            continue
        
        all_genes.update([gene_id])
        probeset = record[probeset_label]
        diff = ExpressionDiff(probeset, fold_change=record[expression_label],
                    prob=record[prob_label],
                    signif=record[sig_label])
        diff.reffile_id = reffile.reffile_id
        diff.sample_a = sample_a
        diff.sample_b = sample_b
        if not successful_commit(session, diff):
            many_probeset += 1
            session.rollback()
            continue
        
        num += 1
    
    session.commit()
    
    run_record.addMessage('add_expression_diff_study',
        LOG_INFO, 'Probeset maps to multiple loci', probeset_to_many)
    run_record.addMessage('add_expression_diff_study',
        LOG_ERROR, 'Probeset in file multiple times', many_probeset)
    run_record.addMessage('add_expression_diff_study',
        LOG_INFO, 'Total probesets added to db', num)
    run_record.addMessage('add_expression_diff_study',
        LOG_INFO, 'Total genes uniquely referenced', len(all_genes))
    
    return run_record

def add_external_genes(session, sample_name, data_path, table, ensembl_release='58',
        ensembl_id_label='ENSEMBL'):
    """adds Expression instances into the database from table

    Arguments:
        - data_path: the reference file path
        - table: the actual expression data table
        - ensembl_release: the Ensembl release
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
    else:
        reffile = reffile[0]
    
    ensembl_ids = table.getRawData(ensembl_id_label)
    genes = session.query(Gene).filter(Gene.ensembl_id.in_(ensembl_ids)).all()
    for gene in genes:
        external = ExternalGene()
        external.gene = gene
        external.reference_file = reffile
        external.sample = sample
        data.append(external)
    
    session.add_all(data)
    session.commit()
