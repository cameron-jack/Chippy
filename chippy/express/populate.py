import datetime
from sqlalchemy import create_engine
from sqlalchemy.exc import IntegrityError

from cogent import LoadTable
from cogent.db.ensembl import HostAccount, Genome

from chippy.express.db import Gene, Transcript, Exon, \
            ExternalGene, Expression, ExpressionDiff, ReferenceFile, Sample, \
            Session, Base, make_session
from chippy.ref.util import chroms
from chippy.parse.r_dump import RDumpToTable

now = datetime.datetime.now()
today = datetime.date(now.year, now.month, now.day)

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
                                ex.Location.Start, ex.Location.End)
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
    samples = []
    for name, description in names_descriptions:
        sample = Sample(name, description)
        samples.append(sample)
    session.add_all(samples)
    session.commit()

def get_transcript_gene_mapping(session, ensembl_release):
    """docstring for get_transcript_gene_mapping"""
    transcript_to_gene = dict(session.query(Transcript.name,
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


def add_expression_study(session, sample_name, data_path, table,
        ensembl_release='58', ensembl_id_label='ENSEMBL',
        expression_label='exp'):
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
    reffile = session.query(ReferenceFile).filter_by(name=data_path).all()
    if len(reffile) == 0:
        reffile = ReferenceFile(data_path, today)
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
    
    sample = session.query(Sample).filter_by(name=sample_name).all()
    assert len(sample) == 1, 'error querying for a sample'
    sample = sample[0]
    
    # order the table in descending order of expression
    table = table.sorted(columns=expression_label, reverse=expression_label)
    rank = 0
    for record in table:
        transcript_ids = record[ensembl_id_label]
        gene_id = single_gene(transcript_to_gene, transcript_ids)
        if gene_id is None:
            continue
        
        expression = Expression(record[expression_label], rank)
        expression.sample = sample
        expression.gene_id = gene_id
        expression.reference_file = reffile
        data.append(expression)
        rank += 1
    
    session.add_all(data)
    session.commit()

def add_expression_diff_study(session, sample_name, data_path, table,
            ref_a_path, ref_b_path,
            ensembl_release='58', ensembl_id_label='ENSEMBL',
            expression_label='exp'):
    """adds Expression instances into the database from table
    
    Arguments:
        - sample_name: the sample name to connect expression instances to
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
    
    ref_a = session.query(ReferenceFile).filter_by(name=ref_a_path).all()
    ref_b = session.query(ReferenceFile).filter_by(name=ref_b_path).all()
    if len(ref_a) == 0 or len(ref_b) == 0:
        raise RuntimeError("Required data not entered yet")
    
    # the reference files
    ref_a = ref_a[0]
    ref_b = ref_b[0]
    data = []
    reffile = session.query(ReferenceFile).filter_by(name=data_path).all()
    if len(reffile) == 0:
        reffile = ReferenceFile(data_path, today)
    else:
        reffile = reffile[0]
    
    # get all gene ID data for the specified Ensembl release
    transcript_to_gene = get_transcript_gene_mapping(session, ensembl_release)
    
    table = table.sorted(columns=expression_label, reverse=expression_label)
    rank = 0
    for record in table:
        transcript_ids = record[ensembl_id_label]
        gene_id = single_gene(transcript_to_gene, transcript_ids)
        if gene_id is None:
            continue
        
        expression_a = session.query(Expression).filter_by(gene_id=gene_id,
                    reference_file=ref_a)
        expression_b = session.query(Expression).filter_by(gene_id=gene_id,
                    reference_file=ref_b)
        
        diff = ExpressionDiff(fold_change=record[expression_label],
                    prob=record[prob_label],
                    signif=record[sig_label])
        diff.reference_file = reffile
        diff.express_A = expression_a
        diff.express_B = expression_b
        data.append(diff)
        rank += 1
    
    session.add_all(data)
    session.commit()
    
