import datetime
from sqlalchemy import create_engine, and_
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm.exc import NoResultFound

from cogent import LoadTable
from cogent.db.ensembl import HostAccount, Genome
from cogent.util.progress_display import display_wrap

from chippy.express.db_schema import Gene, Transcript, Exon, \
            ExternalGene, Expression, ExpressionDiff, ReferenceFile, Sample, \
            Session, Base, make_session
from chippy.ref.util import chroms
from chippy.express.util import single_gene, _one
from chippy.parse.r_dump import RDumpToTable

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
            canonical = gene.CanonicalTranscript.StableId
            for transcript in gene.Transcripts:
                db_transcript = Transcript(transcript.StableId, ensembl_release)
                total_objects += 1
                db_transcript.gene = db_gene
                data.append(db_transcript)
                
                if transcript.StableId == canonical:
                    # add the exons too
                    num_exons = len(transcript.Exons)
                    for ex in transcript.Exons:
                        db_exon = Exon(ex.StableId, ex.Rank,
                                ex.Location.Start, ex.Location.End,
                                ensembl_release)
                        db_exon.transcript = db_transcript
                        db_exon.gene = db_gene
                        data.append(db_exon)
                        total_objects += 1
                
            n += 1
            if n % 100 == 0:
                print 'Genes processed=%s; Db objects created=%d' % (n, total_objects)
                if debug:
                    session.add_all(data)
                    session.commit()
                    return
    
    print 'Writing objects into db'
    session.add_all(data)
    session.commit()


def add_samples(session, names_descriptions):
    """add basic sample type"""
    failed = []
    for name, description in names_descriptions:
        sample = Sample(name, description)
        if not successful_commit(session, sample):
            session.rollback()
            failed.append([name, description])
    
    if len(failed) > 0:
        print 'The following were excluded as they exist in the db'
        print failed

def get_transcript_gene_mapping(session, ensembl_release):
    """docstring for get_transcript_gene_mapping"""
    transcript_to_gene = dict(session.query(Transcript.ensembl_id,
        Transcript.gene_id).filter_by(ensembl_release=ensembl_release).all())
    return transcript_to_gene

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


@display_wrap
def add_expression_study(session, sample_name, data_path, table,
        ensembl_release='58', ensembl_id_label='ENSEMBL',
        expression_label='exp', ui=None):
    """adds Expression instances into the database from table
    
    Arguments:
        - sample_name: the sample name to connect expression instances to
        - data_path: the reference file path
        - table: the actual expression data table
        - ensembl_release: the Ensembl release
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
    
    # check we haven't already added expression data from this file
    records = session.query(Expression).filter_by(reference_file=reffile).all()
    if len(records) > 0:
        print 'Already added this expression file'
        return
    
    # get all gene ID data for the specified Ensembl release
    transcript_to_gene = get_transcript_gene_mapping(session, ensembl_release)
    
    if not successful_commit(session, data):
        pass
    
    data = []
    
    # order the table in descending order of expression
    table = table.sorted(columns=expression_label, reverse=expression_label)
    rank = 0
    failed = []
    probeset_many_loci = 0
    for record in ui.series(table, noun='Adding expression instances'):
        transcript_ids = record[ensembl_id_label]
        gene_id = single_gene(transcript_to_gene, transcript_ids)
        if gene_id is None:
            probeset_many_loci += 1
            continue
        
        expression = Expression(record[expression_label], rank)
        expression.sample = sample
        expression.gene_id = gene_id
        expression.reference_file = reffile
        if not successful_commit(session, expression):
            failed.append(gene_id)
            session.rollback()
            continue
        
        rank += 1
    
    deleted = 0
    if len(failed) > 0:
        failed = session.query(Expression).\
                filter(Expression.gene_id.in_(failed)).all()
        for failure in ui.series(failed,
            noun='Deleting records that had duplicates'):
            session.delete(failure)
            deleted += 1
        
        session.commit()
        # now redo the rank's
        all_expressed = session.query(Expression).filter_by(
                reference_file_id=reffile.reference_file_id).order_by(
                                Expression.rank).all()
        
        for i in ui.series(range(len(all_expressed)),
            noun='Updating ranks for unique expressed records'):
            expressed = all_expressed[i]
            expressed.rank = i
            session.merge(expressed)
        
        session.commit()
    
    loci_many_probesets = len(failed) + deleted
    kept = session.query(Expression).filter_by(
                reference_file_id=reffile.reference_file_id).count()
    print 'Number of probesets excluded because:'
    print '\tProbeset maps to multiple loci = %d' % probeset_many_loci
    print '\tLocus maps to multiple probesets = %d' % loci_many_probesets
    print 'Unique loci with expression: %d' % kept
    

def _one(query):
    """returns result if found, False otherwise"""
    try:
        result = query.one()
    except NoResultFound:
        result = False
    
    return result

@display_wrap
def add_expression_diff_study(session, data_path, table,
            ref_a_path, ref_b_path,
            ensembl_release='58', ensembl_id_label='ENSEMBL',
            expression_label='exp', prob_label='rawp', sig_label='sig',
            ui=None):
    """adds Expression instances into the database from table
    
    Arguments:
        - data_path: the reference file path
        - table: the actual expression data table
        - ref_a_path, ref_b_path: the difference file contains measurements
          between file path a and file path b. This order is CRITICAL as it
          affects the inferences concerning gene expression being up / down
          in a sample.
        - ensembl_release: the Ensembl release
        - ensembl_id_label: label of the column containing Ensembl Stable IDs
        - expression_label: label of the column containing absolute measure
          of expression
        - prob_label: label of the column containing raw probability of
          difference in gene expression between samples A/B
        - sig_label: label of the column classifying probabilities as
          significant after multiple test correction. 1 means up in A relative
          to B, -1 means down in A relative to B, 0 means no difference.
    """
    
    ref_a = _one(session.query(ReferenceFile).filter_by(name=ref_a_path))
    if not ref_a:
        print 'Could not find a record for ref_a %s' % ref_a_path
        print session.query(ReferenceFile).filter_by(name=ref_a_path).all()
    
    ref_b = _one(session.query(ReferenceFile).filter_by(name=ref_b_path))
    if not ref_b:
        print 'Could not find a record for ref_b %s' % ref_b_path
        print session.query(ReferenceFile).filter_by(name=ref_b_path).all()
    
    if not ref_a or not ref_b:
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
    
    data = []
    
    # get all transcript to gene ID data for the specified Ensembl release
    transcript_to_gene = get_transcript_gene_mapping(session, ensembl_release)
    
    table = table.sorted(columns=expression_label, reverse=expression_label)
    rank = 0
    probeset_many_loci = 0
    failed = []
    for record in ui.series(table, noun='Adding expression diffs'):
        transcript_ids = record[ensembl_id_label]
        gene_id = single_gene(transcript_to_gene, transcript_ids)
        if gene_id is None:
            probeset_many_loci += 1
            continue
        
        diff = ExpressionDiff(fold_change=record[expression_label],
                    prob=record[prob_label],
                    signif=record[sig_label])
        diff.reference_file = reffile
        diff.sample_a = sample_a
        diff.sample_b = sample_b
        diff.gene_id = gene_id
        if not successful_commit(session, diff):
            failed.append(gene_id)
            session.rollback()
            continue
        
        rank += 1
    
    deleted = 0
    if len(failed) > 0:
        failed = session.query(ExpressionDiff).\
                filter(ExpressionDiff.gene_id.in_(failed)).all()
        for failure in ui.series(failed,
            noun='Deleting diff records that had duplicates'):
            session.delete(failure)
            deleted += 1
        
        session.commit()
    
    loci_many_probesets = len(failed) + deleted
    kept = session.query(ExpressionDiff).filter_by(
                reference_file_id=reffile.reference_file_id).count()
    print 'Number of probesets excluded because:'
    print '\tProbeset maps to multiple loci = %d' % probeset_many_loci
    print '\tLocus maps to multiple probesets = %d' % loci_many_probesets
    print 'Unique loci with expression diffs: %d' % kept


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
