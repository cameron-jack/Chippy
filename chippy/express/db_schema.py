from __future__ import division
import sys
sys.path.extend(['..'])

from sqlalchemy import (Integer, Float, String, Date, PickleType,
    Boolean, ForeignKey, Column, UniqueConstraint, Table, create_engine)
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import backref, mapper, relationship, sessionmaker

from cogent.util.misc import flatten
from chippy.util.definition import PLUS_STRAND, MINUS_STRAND, NULL_STRAND

__author__ = 'Gavin Huttley, Cameron Jack'
__copyright__ = 'Copyright 2011-2013, Gavin Huttley, Cameron Jack, Anuj Pahwa'
__credits__ = ['Gavin Huttley', 'Cameron Jack']
__license__ = 'GPL'
__maintainer__ = 'Cameron Jack'
__email__ = 'cameron.jack@anu.edu.au'
__status__ = 'Pre-release'
__version__ = '0.1'

###
# IMPORTANT: Coords given by PyCogent and used here are offset 0 [,)
# (Ensembl is offset 1 [,] ) so that when the genome is placed in a numpy
# array, the coords may be used to slice the array directly
###

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
    
    def __str__(self):
        return "%s : %s" % (self.name, self.description)
    

class ReferenceFile(Base):
    """original input source file name"""
    __tablename__ = 'reference_file'
    
    reffile_id = Column(Integer, primary_key=True)
    name = Column(String, unique=True)
    date = Column(Date)
    sample_id = Column(Integer, ForeignKey('sample.sample_id'))
    ref_a_name = Column(String, nullable=True)
    ref_b_name = Column(String, nullable=True)
    
    sample = relationship(Sample,
                backref=backref('reference_files', order_by=reffile_id))
    
    def __init__(self, name, date, ref_a_name=None, ref_b_name=None):
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

class Chroms(Base):
    __tablename__ = 'chroms'
    chroms_id = Column(Integer, primary_key=True)
    species = Column(String, unique=True)
    chromStr = Column(String) # use ',' separated list of chroms

    def __init__(self, species, chromsList):
        self.species = species
        self.chromStr = ','.join(chromsList)

class Gene(Base):
    __tablename__ = 'gene'
    
    gene_id = Column(Integer, primary_key=True)
    
    ensembl_id = Column(String)
    symbol = Column(String)
    biotype = Column(String)
    status = Column(String)
    description = Column(String, nullable=True)
    chrom = Column(String)
    start = Column(Integer)
    end = Column(Integer)
    strand = Column(Integer)
    
    __table_args__ = (UniqueConstraint('ensembl_id', name='unique'), {})
    
    def __init__(self, ensembl_id, symbol, biotype, description, status, chrom, start, end, strand):
        self.ensembl_id = str(ensembl_id)
        self.symbol = symbol
        self.biotype = biotype
        self.description = description
        self.status = status
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self._exon_coords = None
        self._intron_coords = None
        self._tss = None
        self._gene3p = None
        self._probeset_scores = None
        self._rank_stats = None
    
    def __repr__(self):
        return self.ensembl_id, self.chrom, self.start, self.end, self.strand, self.Rank
    
    @property
    def Rank(self):
        if not hasattr(self, '_rank_stats'):
            self._rank_stats = {}
        if not self._rank_stats:
            self._rank_stats = {}

        return self._rank_stats.get('rank', None)

    @Rank.setter
    def Rank(self, value):
        if not hasattr(self, '_rank_stats'):
            self._rank_stats = {}
        
        self._rank_stats['rank'] = value
    
    @property
    def MeanScore(self):
        """returns average expression score of Expression instances"""
        return sum(self.Scores) / len(self.Scores)
    
    @property
    def MaxScore(self):
        """returns maximum score"""
        return max(self._probeset_scores)
    
    @property
    def SumScores(self):
        """sum of expression scores"""
        return sum(self._probeset_scores)
    
    @property
    def Scores(self):
        # getScores
        if not hasattr(self, '_probeset_scores'):
            self._probeset_scores = None
        
        return self._probeset_scores
    
    @Scores.setter
    def Scores(self, value):
        if not hasattr(self, '_probeset_scores'):
            self._probeset_scores = None
        
        self._probeset_scores = flatten(value)
    
    @property
    def IntronCoords(self):
        """returns list of intron coordinates"""
        if not hasattr(self, '_intron_coords'):
            self._intron_coords = None
        
        if self._intron_coords is None:
            exons = self.ExonCoords
            introns = []
            for i in range(1, len(exons)):
                intron = (exons[i-1][1], exons[i][0])
                introns.append(intron)
            self._intron_coords = introns
        
        return self._intron_coords
    
    @property
    def IntronCoordsByRank(self):
        """returns intron coords in their rank order"""
        if self.strand == MINUS_STRAND:
            result = self.IntronCoords[::-1]
        else:
            result = self.IntronCoords
        return result

    @property
    def ExonCoords(self):
        """returns list of exon coordinates"""
        if not hasattr(self, '_exon_coords'):
            self._exon_coords = None
        
        if self._exon_coords is None:
            self._exon_coords = sorted([(e.start, e.end) for e in self.exons])
        return self._exon_coords
    
    @property
    def ExonCoordsByRank(self):
        """returns exon coords in their rank order"""
        if self.strand == MINUS_STRAND:
            result = self.ExonCoords[::-1]
        else:
            result = self.ExonCoords
        return result

    @property
    def Tss(self):
        """the transcription start site"""
        if not hasattr(self, '_tss'):
            self._tss = None
        
        if self._tss is None:
            if self.strand == PLUS_STRAND:
                self._tss = self.start
            else:
                self._tss = self.end
        return self._tss

    @property
    def Gene3p(self):
        """ returns the gene coord at the opposite end to the Tss """
        if not hasattr(self, '_gene3p'):
            self._gene3p = None

        if self._gene3p is None:
            if self.strand == PLUS_STRAND:
                self._gene3p = self.end
            else:
                self._gene3p = self.start
        return self._gene3p

    # All gene related properties and functions appear in 5'->3' order
    # in the code from this point.

    def getTssWindowCoords(self, window_upstream, window_downstream):
        """
            Returns start, finish relative to the gene TSS.
            'end' will be 1 greater than the actual end position
        """
        tss = self.Tss

        if self.strand == PLUS_STRAND:
            start = tss - window_upstream
            end = tss + window_downstream # includes site
        else:
            start = tss - window_downstream
            end = tss + window_upstream # includes site

        assert start < tss < end, 'Start, TSS, End coords' +\
               str(start) + ', ' + str(tss) + ', ' + str(end)

        return start, end

    def getUTRExonWindowCoords(self, window_upstream, window_downstream, no_overlap=True):
        """
            Returns start, finish relative to the 5'UTR/1st Exon boundary.
            Proximity is used to test if this position is within a given distance of
            the TSS. If it is then we exclude these coords.
        """
        # for checking proximity
        tss = self.Tss
        gene3p = self.Gene3p

        if self.strand == PLUS_STRAND:
            exon_start = self.ExonCoordsByRank[0][0] # start of first exon
            start = exon_start - window_upstream
            end = exon_start + window_downstream # includes site
            if (start <= tss or end >= gene3p) and no_overlap:
                return None, None
        else:
            exon_start = self.ExonCoordsByRank[0][1] # start of first exon
            start = exon_start - window_downstream
            end = exon_start + window_upstream # includes site
            if (start <= gene3p or end >= tss) and no_overlap:
                return None, None

        return start, end

    def getIntronExonWindowCoords(self, window_upstream, window_downstream, no_overlap=True):
        """
            Return windows relative to each intron/exon boundary as list.
            We can go by our listing of gene.introns
            Here we are looking at the 'feature' being the 5' end of
            each Intron.
            Introns do not include UTRs.
            Genes with a single exon will return an empty list
            'end' will be 1 greater than the actual end position.
            no_overlap drops IE boundaries are not within the window
            distance of UTR/Exon or Exon/UTR boundaries.
        """

        if len(self.ExonCoordsByRank) < 2:
            return []

        coords = self.ExonCoordsByRank

        UTR5p = self.ExonCoordsByRank[0][0]
        UTR3p = self.ExonCoordsByRank[-1][1]

        all_coords = []
        for c in coords:
            if self.strand == PLUS_STRAND:
                site = c[0]
                start = site - window_upstream
                end = site + window_downstream # includes site
                if (start <= UTR5p or end >= UTR3p) and no_overlap:
                    start = None
                    end = None
            else:
                site = c[1]
                start = site - window_downstream
                end = site + window_upstream # includes site
                if (start >= UTR5p or end <= UTR3p) and no_overlap:
                    start = None
                    end = None
            all_coords.append((start, end))

        return all_coords

    def getExonIntronWindowCoords(self, window_upstream, window_downstream, no_overlap=True):
        """
            Return windows relative to each exon/intron boundary as list.
            We can go by our listing of gene.introns
            Here we are looking at the 'feature' being the 5' end of
            each intron (excluding the first).
            Genes with a single exon will return an empty list
            'end' will be 1 greater than the actual end position.
            no_overlap drops EI boundaries are not within the window
            distance of UTR/Exon or Exon/UTR boundaries.
        """

        if len(self.IntronCoordsByRank) == 0:
            return []

        coords = self.IntronCoordsByRank

        UTR5p = self.ExonCoordsByRank[0][0]
        UTR3p = self.ExonCoordsByRank[-1][1]

        all_coords = []
        for c in coords:
            if self.strand == PLUS_STRAND:
                site = c[0]
                start = site - window_upstream
                end = site + window_downstream # includes site
                if (start <= UTR5p or end >= UTR3p) and no_overlap:
                    start = None
                    end = None
            else:
                site = c[1]
                start = site - window_downstream
                end = site + window_upstream # includes site
                if (start >= UTR5p or end <= UTR3p) and no_overlap:
                    start = None
                    end = None
            all_coords.append((start, end))

        return all_coords

    def getExonUTRWindowCoords(self, window_upstream, window_downstream, no_overlap=True):
        """
            Returns start, finish relative to the last Exon/3'UTR boundary.
            Proximity is used to test if this position is within a given distance of
            the 3' end of the gene. If it is then we exclude these coords.
        """
        # for checking proximity
        tss = self.Tss
        gene3p = self.Gene3p

        if self.strand == PLUS_STRAND:
            exon_end = self.ExonCoordsByRank[-1][1] # start of 3' UTR
            start = exon_end - window_upstream
            end = exon_end + window_downstream # includes site
            if (start <= tss or end >= gene3p) and no_overlap:
                return None, None
        else:
            exon_end = self.ExonCoordsByRank[-1][0] # start of 3' UTR
            start = exon_end - window_downstream
            end = exon_end + window_upstream # includes site
            if (start <= gene3p or end >= tss) and no_overlap:
                return None, None

        return start, end

    def getGene3PrimeWindowCoords(self, window_upstream, window_downstream):
        """
            Return opposite end of the gene to the TSS
        """
        end_site = self.Gene3p

        if self.strand == PLUS_STRAND:
            start = end_site - window_upstream
            end = end_site + window_downstream # includes site
        else:
            start = end_site - window_downstream
            end = end_site + window_upstream # includes site

        assert start < end_site < end, 'Start, Gene-3p, End coords' +\
                          str(start) + ', ' + str(end_site) + ', ' + str(end)

        return start, end

class Exon(Base):
    __tablename__ = 'exon'
    
    exon_id = Column(Integer, primary_key=True)
    
    ensembl_id = Column(String)
    rank = Column(Integer)
    start = Column(Integer)
    end = Column(Integer)
    
    gene_id = Column(Integer, ForeignKey('gene.gene_id'))
    gene = relationship(Gene,
                backref=backref('exons', order_by=exon_id))
    
    __table_args__ = (UniqueConstraint('gene_id', 'rank', 'ensembl_id',
                        name='unique'), {})
    
    
    def __init__(self, ensembl_id, rank, start, end):
        super(Exon, self).__init__()
        self.ensembl_id = ensembl_id
        self.start = start
        self.end = end
        self.rank = rank
    
    def __repr__(self):
        return "Exon(gene=%s, start=%s, rank=%s, strand=%s)" % (
            self.gene.ensembl_id, self.start, self.rank,
            self.gene.strand)

class Expression(Base):
    __tablename__ = 'expression'
    
    expression_id = Column(Integer, primary_key=True)
    
    probesets = Column(PickleType)
    scores = Column(PickleType)
    
    sample_id = Column(Integer, ForeignKey('sample.sample_id'))
    gene_id = Column(Integer, ForeignKey('gene.gene_id'))
    reffile_id = Column(Integer,
            ForeignKey('reference_file.reffile_id'))
    
    sample = relationship(Sample,
                backref=backref('expression', order_by=expression_id))
    gene = relationship(Gene, backref=backref('expression'))
    
    __table_args__ = (UniqueConstraint('gene_id', 'sample_id', 'reffile_id',
                        name='unique'), {})
    
    def __init__(self, probesets, scores):
        super(Expression, self).__init__()
        self.probesets = probesets
        self.scores = scores
    
    def __repr__(self):
        return 'Expression(probesets=%s, sample=%s)' % (
                self.probesets, self.sample.name)

class ExpressionDiff(Base):
    __tablename__ = 'expression_diff'
    
    expression_diff_id = Column(Integer, primary_key=True)
    
    probesets = Column(PickleType)
    fold_changes = Column(PickleType)
    probability = Column(Float)
    multitest_signif = Column(Integer)
    
    gene_id = Column(Integer, ForeignKey('gene.gene_id'))
    sample_id = Column(Integer, ForeignKey('sample.sample_id'))
    reffile_id = Column(Integer, ForeignKey('reference_file.reffile_id'))
    
    sample = relationship(Sample,
                backref=backref('expression_diff', order_by=expression_diff_id))
    
    gene = relationship(Gene,
            backref=backref('expression_diff', order_by=expression_diff_id))
    
    reference_file = relationship(ReferenceFile,
            backref=backref('expression_diff', order_by=expression_diff_id))
    
    __table_args__ = (UniqueConstraint('gene_id', 'reffile_id',
                        name='unique'), {})
    
    def __init__(self, probesets, fold_changes, prob, signif):
        super(ExpressionDiff, self).__init__()
        self.probesets = probesets
        self.fold_changes = fold_changes
        self.probability = prob
        self.multitest_signif = signif
    
    def __repr__(self):
        return 'ExpressionDiff(probesets=%s, P=%s, Signif=%s)' %\
            (self.probesets, self.probability, self.multitest_signif)


class TargetGene(Base):
    """a gene of interest without expression info"""
    __tablename__ = 'target_gene'
    target_gene_id = Column(Integer, primary_key=True)
    
    rank = Column(Integer)
    
    sample_id = Column(Integer, ForeignKey('sample.sample_id'))
    gene_id = Column(Integer, ForeignKey('gene.gene_id'))
    reffile_id = Column(Integer,
            ForeignKey('reference_file.reffile_id'))
    
    sample = relationship(Sample,
                backref=backref('target_gene', order_by=target_gene_id))
    gene = relationship(Gene,
                backref=backref('target_gene', order_by=target_gene_id))
    reference_file = relationship(ReferenceFile,
                backref=backref('target_gene', order_by=target_gene_id))
    
    __table_args__ = (UniqueConstraint('gene_id', 'reffile_id',
                    name='unique'), {})



def _make_session(db_name, reset=False):
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