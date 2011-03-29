from sqlalchemy import (Integer, Float, String, Date,
    Boolean, ForeignKey, Column, UniqueConstraint, Table, create_engine)
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import backref, mapper, relationship, sessionmaker

Session = sessionmaker()

Base = declarative_base()

class Sample(Base):
    __tablename__ = "sample"
    
    sample_id = Column(Integer, primary_key=True)
    
    name = Column(String, unique=True)
    description = Column(String)
    
    def __init__(self, name, description):
        self.name = name
        self.description = description
    
    def __repr__(self):
        return "Sample('%s', '%s')" % (self.name, self.description)
    

class ReferenceFile(Base):
    """original input source file name"""
    __tablename__ = 'reference_file'
    
    reference_file_id = Column(Integer, primary_key=True)
    name = Column(String, unique=True)
    date = Column(Date)
    sample_id = Column(Integer, ForeignKey('sample.sample_id'))
    ref_a_name = Column(String, nullable=True)
    ref_b_name = Column(String, nullable=True)
    
    sample = relationship(Sample,
                backref=backref('reference_files', order_by=reference_file_id))
    
    def __init__(self, name,date, ref_a_name=None, ref_b_name=None):
        super(ReferenceFile, self).__init__()
        self.name = name
        self.date = date
        if ref_a_name or ref_b_name:
            assert ref_a_name and ref_b_name, 'Need 2 reference file names'
        
        self.ref_a_name = ref_a_name
        self.ref_b_name = ref_b_name
    
    def __repr__(self):
        depends = ''
        if self.ref_a_name:
            depends = ', depends=%s' % str([self.ref_a_name, self.ref_b_name])
        return "ReferenceFile('%s'%s)" % (self.name, depends)
    

class Gene(Base):
    __tablename__ = 'gene'
    
    gene_id = Column(Integer, primary_key=True)
    
    ensembl_id = Column(String)
    ensembl_release = Column(String)
    symbol = Column(String)
    biotype = Column(String)
    status = Column(String)
    description = Column(String, nullable=True)
    coord_name = Column(String)
    start = Column(Integer)
    end = Column(Integer)
    strand = Column(Integer)
    
    __table_args__ = (UniqueConstraint('ensembl_id', 'ensembl_release',
                name='unique'), {})
    
    def __init__(self, ensembl_release, ensembl_id, symbol, biotype, description, status, coord_name, start, end, strand):
        self.ensembl_release = ensembl_release
        self.ensembl_id = ensembl_id
        
        self.symbol = symbol
        self.biotype = biotype
        self.description = description
        self.status = status
        self.coord_name = coord_name
        self.start = start
        self.end = end
        self.strand = strand
        self._exon_coords = None
        self._intron_coords = None
        self._tss = None
    
    def __repr__(self):
        return "Gene(ensembl_id='%s', coord_name='%s', start=%s, strand=%s)" \
                % (self.ensembl_id, self.coord_name, self.start, self.strand)
    
    @property
    def IntronCoords(self):
        """returns list of intron coordinates"""
        if self._intron_coords is None:
            exons = self.ExonCoords
            introns = []
            for i in range(1, len(exons)):
                intron = (exons[i-1][1], exons[i][0])
                introns.append(intron)
            self._intron_coords = introns
        
        return self._intron_coords
    
    @property
    def ExonCoords(self):
        """returns list of exon coordinates"""
        if self._exon_coords is None:
            self._exon_coords = [(e.start, e.end) for e in self.exons]
        return self._exon_coords
    
    @property
    def Tss(self):
        """the transcription start site"""
        if self._tss is None:
            if self.strand == 1:
                self._tss = self.start
            else:
                self._tss = self.end
        return self._tss
    
    def getTssCentredCoords(self, size):
        """returns coords centred on the gene TSS"""
        tss = self.Tss
        return sorted([tss-size, tss+size])
    
    def getUpstreamCoords(self, size):
        """returns coords ending at the TSS"""
        tss = self.Tss
        if self.strand == 1:
            start, end = tss-size, tss
        else:
            start, end = tss, tss+size
        return start, end
    

class Transcript(Base):
    """map Ensembl transcripts to genes"""
    __tablename__ = 'transcript'
    
    transcript_id = Column(Integer, primary_key=True)
    
    ensembl_id = Column(String, unique=True)
    ensembl_release = Column(String)
    
    gene_id = Column(Integer, ForeignKey('gene.gene_id'))
    gene = relationship(Gene,
                backref=backref('transcripts', order_by=transcript_id))
    
    __table_args__ = (UniqueConstraint('ensembl_id', 'ensembl_release',
                name='unique'), {})
    
    def __init__(self, ensembl_id, ensembl_release):
        super(Transcript, self).__init__()
        self.ensembl_id = ensembl_id
        self.ensembl_release = ensembl_release
        
    def __repr__(self):
        return "Transcript(%s, %s)" % (self.ensembl_id, self.gene.ensembl_id)
    

class Exon(Base):
    __tablename__ = 'exon'
    
    exon_id = Column(Integer, primary_key=True)
    
    ensembl_id = Column(String)
    rank = Column(Integer)
    start = Column(Integer)
    end = Column(Integer)
    ensembl_release = Column(String)
    
    transcript_id = Column(Integer, ForeignKey('transcript.transcript_id'))
    transcript = relationship(Transcript,
                backref=backref('exons', order_by=exon_id))
    
    gene_id = Column(Integer, ForeignKey('gene.gene_id'))
    gene = relationship(Gene,
                backref=backref('exons', order_by=rank))
    
    __table_args__ = (UniqueConstraint('transcript_id', 'rank', 'ensembl_id',
                        name='unique'), {})
    
    
    def __init__(self, ensembl_id, rank, start, end, ensembl_release):
        super(Exon, self).__init__()
        self.ensembl_id = ensembl_id
        self.start = start
        self.end = end
        self.rank = rank
        self.ensembl_release = ensembl_release
    
    def __repr__(self):
        return "Exon(gene=%s, start=%s, rank=%s, strand=%s)" % (
            self.transcript.gene.ensembl_id, self.start, self.rank,
            self.gene.strand)
        
    

class Expression(Base):
    __tablename__ = 'expression'
    
    expression_id = Column(Integer, primary_key=True)
    
    expression_score = Column(Float)
    rank = Column(Integer)
    
    sample_id = Column(Integer, ForeignKey('sample.sample_id'))
    gene_id = Column(Integer, ForeignKey('gene.gene_id'))
    reference_file_id = Column(Integer,
            ForeignKey('reference_file.reference_file_id'))
    
    sample = relationship(Sample,
                backref=backref('expression', order_by=expression_id))
    gene = relationship(Gene,
                backref=backref('expression', order_by=expression_id))
    reference_file = relationship(ReferenceFile,
                backref=backref('expression', order_by=expression_id))
    
    __table_args__ = (UniqueConstraint('gene_id', 'reference_file_id',
                        name='unique'), {})
    
    def __init__(self, expression_score, rank):
        super(Expression, self).__init__()
        self.expression_score = expression_score
        self.rank = rank
    
    def __repr__(self):
        return 'Expression(ensembl_id=%s, sample=%s, score=%s, rank=%s)' % (
                self.gene.ensembl_id, self.sample.name, self.expression_score,
                self.rank)
        
    


class ExpressionDiff(Base):
    __tablename__ = 'expression_diff'
    
    expression_diff_id = Column(Integer, primary_key=True)
    
    fold_change = Column(Float)
    probability = Column(Float)
    multitest_signif = Column(Integer)
    
    gene_id = Column(Integer, ForeignKey('gene.gene_id'))
    sample_a_id = Column(Integer, ForeignKey('sample.sample_id'))
    sample_b_id = Column(Integer, ForeignKey('sample.sample_id'))
    reference_file_id = Column(Integer,
            ForeignKey('reference_file.reference_file_id'))
    
    gene = relationship(Gene,
                backref=backref('expression_diffs', order_by=expression_diff_id))
    sample_a = relationship(Sample,
            primaryjoin = sample_a_id == Sample.sample_id)
    sample_b = relationship(Sample,
            primaryjoin = sample_b_id == Sample.sample_id)
    
    reference_file = relationship(ReferenceFile,
            backref=backref('expression_diffs', order_by=expression_diff_id))
    
    __table_args__ = (UniqueConstraint('gene_id', 'reference_file_id',
                        name='unique'), {})
    
    
    def __init__(self, fold_change, prob, signif):
        super(ExpressionDiff, self).__init__()
        self.fold_change = fold_change
        self.probability = prob
        self.multitest_signif = signif
    
    def __repr__(self):
        return 'ExpressionDiff(ensembl_id=%s, A=%s, B=%s, P=%s, Signif=%s)' %\
            (self.gene.ensembl_id, self.sample_a.name,
            self.sample_b.name, self.probability, self.multitest_signif)


class ExternalGene(Base):
    """a gene of interest identified by an external source"""
    __tablename__ = 'external_gene'
    external_gene_id = Column(Integer, primary_key=True)
    
    rank = Column(Integer)
    
    sample_id = Column(Integer, ForeignKey('sample.sample_id'))
    gene_id = Column(Integer, ForeignKey('gene.gene_id'))
    reference_file_id = Column(Integer,
            ForeignKey('reference_file.reference_file_id'))
    
    sample = relationship(Sample,
                backref=backref('external_genes', order_by=external_gene_id))
    gene = relationship(Gene,
                backref=backref('external_gene', order_by=external_gene_id))
    reference_file = relationship(ReferenceFile,
                backref=backref('external_genes', order_by=external_gene_id))
    
    __table_args__ = (UniqueConstraint('gene_id', 'reference_file_id',
                    name='unique'), {})
    

def make_session(db_name, reset=False):
    """returns a db session given the db_name"""
    engine = create_engine(db_name)
    metadata = Base.metadata
    if reset:
        metadata.drop_all(engine)
        pass
    
    metadata.create_all(engine)
    Session.configure(bind=engine)
    session = Session()
    return session
