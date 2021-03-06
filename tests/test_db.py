"""test design of the express.db module"""
import sys
sys.path.append('..')

import warnings
warnings.filterwarnings('ignore', 'Not using MPI as mpi4py not found')

from cogent.util.unit_test import TestCase, main

import datetime
from sqlalchemy.exc import IntegrityError

from chippy.express.db_schema import Gene, Exon, TargetGene, \
        Expression, ExpressionDiff, ReferenceFile, Sample

from chippy.express.db_query import get_sample_counts, get_gene_counts, \
        get_expression_counts, get_diff_counts, get_targetgene_counts, \
        get_reffile_counts, get_exon_counts, \
        get_sample_entries, get_gene_entries, get_expression_entries, \
        get_diff_entries, get_targetgene_entries, get_reffile_entries, \
        get_exon_entries, \
        get_genes_by_ranked_expr, get_genes_by_ranked_diff, get_chroms, \
        get_species, make_session, drop_sample_records

from chippy.express.db_populate import add_expression_diff_study, \
        add_sample, add_chroms
from chippy.parse.expr_data import gene_expr_to_table, gene_expr_diff_to_table

__author__ = 'Gavin Huttley, Cameron Jack'
__copyright__ = 'Copyright 2011-2013, Gavin Huttley, Cameron Jack, Anuj Pahwa'
__credits__ = ['Gavin Huttley', 'Cameron Jack']
__license__ = 'GPL'
__maintainer__ = 'Cameron Jack'
__email__ = 'cameron.jack@anu.edu.au'
__status__ = 'Pre-release'
__version__ = '0.2'

now = datetime.datetime.now()
today = datetime.date(now.year, now.month, now.day)

def add_all_gene_exons(session, genes):
    data = []
    for record in genes:
        gene_data = record['gene']
        exons_data = record['exons']
        gene = Gene(**gene_data)
        data.append(gene)
        for e_data in exons_data:
            exon = Exon(**e_data)
            exon.gene = gene
            data.append(exon)
    session.add_all(data)
    session.commit()

class TestDbBase(TestCase):
    def setUp(self):
        self.session = make_session(":memory:")

class TestRefFiles(TestDbBase):
    a = 'reffile-a.txt'
    b = 'reffile-b.txt'
    d = 'reffile-depends.txt'
    
    def test_depends(self):
        """a reference file with dependencies should correctly link"""
        reffile_a = ReferenceFile(self.a, today)
        reffile_b = ReferenceFile(self.b, today)
        reffile_d = ReferenceFile(self.d, today, ref_a_name=self.a,
                                ref_b_name=self.b)
        self.assertEqual(str(reffile_d),
        "ReferenceFile('reffile-depends.txt', depends=['reffile-a.txt', 'reffile-b.txt'])")
    
    def test_sample_association(self):
        """correctly associate a reference file with a sample"""
        sample = Sample('A', 'a sample')
        reffile_a = ReferenceFile(self.a, today)
        reffile_a.sample = sample
        self.session.add_all([reffile_a, sample])
        self.session.commit()
        reffiles = get_reffile_entries(self.session)
        self.assertEqual(reffiles[0].sample.name, 'A')
        num_reffiles = get_reffile_counts(self.session)
        self.assertEqual(num_reffiles, 1)

class TestChrom(TestDbBase):
    """ correctly set & get chromosomes for a species """
    chromlist = ['1','2','3','4','X','Y']
    species = 'Artificial'

    def test_add_and_get_species_chroms(self):
        self.assertTrue(add_chroms(self.session, self.species, self.chromlist))

        stored_chroms = get_chroms(self.session)
        self.assertTrue(set(self.chromlist)==set(stored_chroms))
        stored_species = get_species(self.session)
        self.assertTrue(stored_species == self.species)

class TestGene(TestDbBase):
    """test gene properties"""
    plus_coords_one_exons = dict(gene=dict(ensembl_id='PLUS-1',
        symbol='agene', biotype='protein_coding', status='fake',
        description='a fake gene',
        chrom='1', start=1000, end=2000, strand=1),
        exons=[dict(ensembl_id='exon-1', rank=1, start=1050, end=1950)]
        )
    
    plus_coords_many_exons = dict(gene=dict(ensembl_id='PLUS-3',
        symbol='agene', biotype='protein_coding', status='fake',
        description='a fake gene',
        chrom='2', start=1000, end=2000, strand=1),
        exons=[dict(ensembl_id='exon-1', rank=1, start=1050, end=1400),
               dict(ensembl_id='exon-2', rank=2, start=1600, end=1700),
               dict(ensembl_id='exon-3', rank=3, start=1800, end=1900)]
        )
    
    # 
    minus_coords_one_exons = dict(gene=dict(ensembl_id='MINUS-1',
        symbol='agene', biotype='protein_coding', status='fake',
        description='a fake gene',
        chrom='2', start=1000, end=2000, strand=-1),
        exons=[dict(ensembl_id='exon-1', rank=1, start=1050, end=1950)]
        )
    
    minus_coords_many_exons = dict(gene=dict(ensembl_id='MINUS-3',
        symbol='agene', biotype='protein_coding', status='fake',
        description='a fake gene',
        chrom='3', start=1000, end=2000, strand=-1),
        exons=[dict(ensembl_id='exon-3', rank=3, start=1050, end=1400),
               dict(ensembl_id='exon-2', rank=2, start=1600, end=1700),
               dict(ensembl_id='exon-1', rank=1, start=1800, end=1900)]
        )
    
    genes = [plus_coords_one_exons, plus_coords_many_exons,
            minus_coords_one_exons, minus_coords_many_exons]
    
    def test_add_genes(self):
        """exercise adding a gene"""
        data = [Gene(**self.plus_coords_many_exons['gene']),
                Gene(**self.plus_coords_one_exons['gene']),
                Gene(**self.minus_coords_many_exons['gene']),
                Gene(**self.minus_coords_one_exons['gene'])]
        
        self.session.add_all(data)
        self.session.commit()
    
    def test_unique_constraint_gene(self):
        """adding same gene/ensembl release should raise IntegrityError"""
        data = [Gene(**self.plus_coords_many_exons['gene']),
                Gene(**self.plus_coords_many_exons['gene']),
                Gene(**self.plus_coords_one_exons['gene'])]
        
        self.session.add_all(data)
        self.assertRaises(IntegrityError, self.session.commit)
    
    def test_unique_constraint_exon(self):
        """adding same exon/rank for a gene should raise IntegrityError"""
        data = []
        gene = Gene(**self.plus_coords_many_exons['gene'])
        data.append(gene)
        for e_data in self.plus_coords_many_exons['exons']:
            exon = Exon(**e_data)
            exon.gene = gene
            data.append(exon)
            
        for e_data in self.plus_coords_many_exons['exons']:
            exon = Exon(**e_data)
            exon.gene = gene
            data.append(exon)
        
        self.session.add_all(data)
        self.assertRaises(IntegrityError, self.session.commit)

    def test_get_exon_entries(self):
        """ correctly return exons """
        add_all_gene_exons(self.session, self.genes)

        exons = get_exon_entries(self.session)
        starts = [exon.start for exon in exons]
        expected = [1050, 1050, 1600, 1800, 1050, 1050, 1600, 1800]
        self.assertEqual(starts, expected)

        # by gene
        exons = get_exon_entries(self.session, gene_ensembl='MINUS-3')
        starts = [exon.start for exon in exons]
        expected = [1800, 1600, 1050]
        self.assertEqual(starts, expected)

        # by gene that doesn't exist
        exons = get_exon_entries(self.session, gene_ensembl='Does_not_exist')
        self.assertEqual(exons, [])

    def test_get_exon_counts(self):
        """ correctly return exon counts """
        add_all_gene_exons(self.session, self.genes)

        num_exons = get_exon_counts(self.session)
        self.assertEqual(num_exons, 8)

        # by gene
        num_exons = get_exon_counts(self.session, gene_ensembl='MINUS-3')
        self.assertEqual(num_exons, 3)

        # by gene that doesn't exist
        num_exons = get_exon_counts(self.session, gene_ensembl='Does_not_exist')
        self.assertEqual(num_exons, 0)
    
    def test_get_gene_exon_coords(self):
        """Gene instances correctly derive coords for their exons"""
        add_all_gene_exons(self.session, self.genes)
        genes = self.session.query(Gene).all()
        
        expect = {'PLUS-1': [(1050, 1950)],
            'PLUS-3': [(1050, 1400),(1600, 1700),(1800, 1900)],
            'MINUS-1': [(1050, 1950)],
            'MINUS-3': [(1050, 1400),(1600, 1700),(1800, 1900)],}
        for gene in genes:
            self.assertEqual(gene.ExonCoords, expect[gene.ensembl_id])
    
    def test_get_gene_exon_coords_by_rank(self):
        """return exon coordinates in correct exon.rank order"""
        add_all_gene_exons(self.session, self.genes)
        genes = self.session.query(Gene).all()
        
        expect = {'PLUS-1': [(1050, 1950)],
            'PLUS-3': [(1050, 1400),(1600, 1700),(1800, 1900)],
            'MINUS-1': [(1050, 1950)],
            'MINUS-3': [(1800, 1900),(1600, 1700),(1050, 1400)],}
        
        for gene in genes:
            self.assertEqual(gene.ExonCoordsByRank, expect[gene.ensembl_id])
    
    def test_get_gene_intron_coords(self):
        """Gene instances correctly derive coords for their intron"""
        add_all_gene_exons(self.session, self.genes)
        genes = self.session.query(Gene).all()
        expect = {'PLUS-1': [],
            'PLUS-3': [(1400, 1600),(1700, 1800)],
            'MINUS-1': [],
            'MINUS-3': [(1400, 1600),(1700, 1800)],}
        for gene in genes:
            self.assertEqual(gene.IntronCoords, expect[gene.ensembl_id])
    
    def test_get_gene_intron_coords_by_rank(self):
        """get intron coords by exon.rank"""
        add_all_gene_exons(self.session, self.genes)
        genes = self.session.query(Gene).all()
        expect = {'PLUS-1': [],
            'PLUS-3': [(1400, 1600),(1700, 1800)],
            'MINUS-1': [],
            'MINUS-3': [(1700, 1800),(1400, 1600)],}
        for gene in genes:
            self.assertEqual(gene.IntronCoordsByRank, expect[gene.ensembl_id])
    
    def test_gene_tss(self):
        """return correct TSS coordinate"""
        add_all_gene_exons(self.session, self.genes)
        genes = self.session.query(Gene).all()
        expect = {'PLUS-1': 1000,
            'PLUS-3': 1000,
            'MINUS-1': 2000,
            'MINUS-3': 2000}
        for gene in genes:
            self.assertEqual(gene.Tss, expect[gene.ensembl_id])

    def test_gene_3_prime(self):
        """ return opposite gene end coord to TSS """
        add_all_gene_exons(self.session, self.genes)
        genes = self.session.query(Gene).all()
        expect = {'PLUS-1': 2000,
                  'PLUS-3': 2000,
                  'MINUS-1': 1000,
                  'MINUS-3': 1000}
        for gene in genes:
            self.assertEqual(gene.Gene3p, expect[gene.ensembl_id])

    def test_getTssWindowCoords(self):
        """
            Return correct coords up and down-stream from TSS.
        """
        add_all_gene_exons(self.session, self.genes)
        genes = self.session.query(Gene).all()

        # symmetrical tests first
        upstream = 200
        downstream = 200
        expect = {'PLUS-1': (1000 - upstream, 1000 + downstream),
                  'PLUS-3': (1000 - upstream, 1000 + downstream),
                  'MINUS-1': (2000 - downstream, 2000 + upstream),
                  'MINUS-3': (2000 - downstream, 2000 + upstream)}
        for gene in genes:
            region = gene.getTssWindowCoords(upstream, downstream)
            self.assertEqual(region, expect[gene.ensembl_id])

        # asymmetric
        upstream = 200
        downstream = 100
        expect = {'PLUS-1': (1000 - upstream, 1000 + downstream),
                  'PLUS-3': (1000 - upstream, 1000 + downstream),
                  'MINUS-1': (2000 - downstream, 2000 + upstream),
                  'MINUS-3': (2000 - downstream, 2000 + upstream)}
        for gene in genes:
            region = gene.getTssWindowCoords(upstream, downstream)
            self.assertEqual(region, expect[gene.ensembl_id])

    def test_getUTRExonWindowCoords(self):
        """
            Returns start, finish relative to the 5'UTR/1st Exon boundary.
            Proximity is used to test if this position is within a given
            distance of the TSS. If it is then we exclude these coords.
        """
        add_all_gene_exons(self.session, self.genes)
        genes = self.session.query(Gene).all()

        # symmetric window
        upstream = 200
        downstream = 200
        expect = {'PLUS-1': (1050 - upstream, 1050 + downstream),
                  'PLUS-3': (1050 - upstream, 1050 + downstream),
                  'MINUS-1': (1950 - downstream, 1950 + upstream),
                  'MINUS-3': (1900 - downstream, 1900 + upstream)}
        for gene in genes:
            region = gene.getUTRExonWindowCoords(upstream, downstream, no_overlap=False)
            self.assertEqual(region, expect[gene.ensembl_id])

        # asymmetric window
        upstream = 200
        downstream = 100
        expect = {'PLUS-1': (1050 - upstream, 1050 + downstream),
                  'PLUS-3': (1050 - upstream, 1050 + downstream),
                  'MINUS-1': (1950 - downstream, 1950 + upstream),
                  'MINUS-3': (1900 - downstream, 1900 + upstream)}
        for gene in genes:
            region = gene.getUTRExonWindowCoords(upstream, downstream, no_overlap=False)
            self.assertEqual(region, expect[gene.ensembl_id])

        # now check for a clash with the TSS or Gene3p site
        upstream = 50
        downstream = 100
        expect = {'PLUS-1': (None, None),
                  'PLUS-3': (None, None),
                  'MINUS-1': (None, None),
                  'MINUS-3': (1900 - downstream, 1900 + upstream)}
        for gene in genes:
            region = gene.getUTRExonWindowCoords(upstream, downstream, no_overlap=True)
            self.assertEqual(region, expect[gene.ensembl_id])


    def test_getExonUTRWindowCoords(self):
        """
            Returns start, finish relative to the last Exon/3'UTR boundary.
            Proximity is used to test if this position is within a given
            distance of the 3' end of the gene. If it is then we exclude
            these coords.
        """
        add_all_gene_exons(self.session, self.genes)
        genes = self.session.query(Gene).all()

        # symmetric window
        upstream = 200
        downstream = 200
        expect = {'PLUS-1': (1950 - upstream, 1950 + downstream),
                  'PLUS-3': (1900 - upstream, 1900 + downstream),
                  'MINUS-1': (1050 - downstream, 1050 + upstream),
                  'MINUS-3': (1050 - downstream, 1050 + upstream)}
        for gene in genes:
            region = gene.getExonUTRWindowCoords(upstream, downstream, no_overlap=False)
            self.assertEqual(region, expect[gene.ensembl_id])

        # asymmetric window
        upstream = 200
        downstream = 100
        expect = {'PLUS-1': (1950 - upstream, 1950 + downstream),
                  'PLUS-3': (1900 - upstream, 1900 + downstream),
                  'MINUS-1': (1050 - downstream, 1050 + upstream),
                  'MINUS-3': (1050 - downstream, 1050 + upstream)}
        for gene in genes:
            region = gene.getExonUTRWindowCoords(upstream, downstream, no_overlap=False)
            self.assertEqual(region, expect[gene.ensembl_id])

        # now check for a clash with the TSS or Gene3p site
        upstream = 100
        downstream = 50
        expect = {'PLUS-1': (None, None),
                  'PLUS-3': (1900 - upstream, 1900 + downstream),
                  'MINUS-1': (None, None),
                  'MINUS-3': (None, None)}
        for gene in genes:
            region = gene.getExonUTRWindowCoords(upstream, downstream, no_overlap=True)
            self.assertEqual(region, expect[gene.ensembl_id])


    def test_getIntronExonWindowCoords(self):
        """
            Return correct coords up and down-stream from Exon5p start.
            Single exon genes will return no boundaries.
            The first exon will be ignored as it borders the UTR
        """

        add_all_gene_exons(self.session, self.genes)
        plus1 = self.session.query(Gene).filter_by(ensembl_id='PLUS-1').one()
        plus3 = self.session.query(Gene).filter_by(ensembl_id='PLUS-3').one()
        minus1 = self.session.query(Gene).filter_by(ensembl_id='MINUS-1').one()
        minus3 = self.session.query(Gene).filter_by(ensembl_id='MINUS-3').one()

        # symmetric window first
        upstream = 100
        downstream = 100

        expect = {'PLUS-1': [], 'PLUS-3': [(c-upstream, c+downstream)
                                                for c in [1050, 1600, 1800]],
                 'MINUS-1': [], 'MINUS-3': [(c-downstream, c+upstream)
                                                for c in [1900, 1700, 1400]]}

        for gene in (plus1, plus3, minus1, minus3):
            got = gene.getIntronExonWindowCoords(upstream, downstream, no_overlap=False)
            self.assertEqual(got, expect[gene.ensembl_id])

        # asymmetric window: 200 upstream, 100 downstream
        upstream = 200
        downstream = 100
        expect = {'PLUS-1': [], 'PLUS-3': [(c-upstream, c+downstream)
                                                for c in [1050, 1600, 1800]],
                  'MINUS-1': [], 'MINUS-3': [(c-downstream, c+upstream)
                                                for c in [1900, 1700, 1400]]}
        for gene in (plus1, plus3, minus1, minus3):
            got = gene.getIntronExonWindowCoords(upstream, downstream, no_overlap=False)
            self.assertEqual(got, expect[gene.ensembl_id])

        # proximity check - wipes out all intron starts close to UTR

        expect = {'PLUS-1': [], 'PLUS-3': [
                    (None, None),
                    (1600-upstream, 1600+downstream),
                    (None, None)
                ], 'MINUS-1': [], 'MINUS-3': [
                (None, None),
                (1700-downstream, 1700+upstream),
                (1400-downstream, 1400+upstream)
                ]}

        for gene in (plus1, plus3, minus1, minus3):
            got = gene.getIntronExonWindowCoords(upstream, downstream, no_overlap=True)
            self.assertEqual(got, expect[gene.ensembl_id])


    def test_getExonIntronWindowCoords(self):
        """
            Return correct coords up and down-stream from Intron5p start
            Single exon genes will return no boundaries.
            The last exon will be ignored as it borders UTR.
        """
        add_all_gene_exons(self.session, self.genes)
        plus1 = self.session.query(Gene).filter_by(ensembl_id='PLUS-1').one()
        plus3 = self.session.query(Gene).filter_by(ensembl_id='PLUS-3').one()
        minus1 = self.session.query(Gene).filter_by(ensembl_id='MINUS-1').one()
        minus3 = self.session.query(Gene).filter_by(ensembl_id='MINUS-3').one()

        # symmetric window first
        upstream = 100
        downstream = 100

        expect = {'PLUS-1': [], 'PLUS-3': [(c-upstream, c+downstream)
                                                for c in [1400, 1700]],
                  'MINUS-1': [], 'MINUS-3': [(c-downstream, c+upstream)
                                                for c in [1800, 1600]]}

        for gene in (plus1, plus3, minus1, minus3):
            got = gene.getExonIntronWindowCoords(upstream, downstream, no_overlap=False)
            self.assertEqual(got, expect[gene.ensembl_id])

        # asymmetric window: 200 upstream, 100 downstream
        upstream = 200
        downstream = 100
        expect = {'PLUS-1': [], 'PLUS-3': [(c-upstream, c+downstream)
                                            for c in [1400, 1700]],
                  'MINUS-1': [], 'MINUS-3': [(c-downstream, c+upstream)
                                            for c in [1800, 1600]]}

        for gene in (plus1, plus3, minus1, minus3):
            got = gene.getExonIntronWindowCoords(upstream, downstream, no_overlap=False)
            self.assertEqual(got, expect[gene.ensembl_id])

        # proximity check - wipes out all intron starts close to UTR

        expect = {'PLUS-1': [], 'PLUS-3': [(c-upstream, c+downstream)
                                            for c in [1400, 1700]],
                  'MINUS-1': [], 'MINUS-3': [(c-downstream, c+upstream)
                                            for c in [1800, 1600]]}

        for gene in (plus1, plus3, minus1, minus3):
            got = gene.getExonIntronWindowCoords(upstream, downstream, no_overlap=True)
            self.assertEqual(got, expect[gene.ensembl_id])

    def test_getGene3PrimeWindowCoords(self):
        """
            Return correct coords up and down-stream from gene 3prime location
        """
        add_all_gene_exons(self.session, self.genes)
        genes = self.session.query(Gene).all()

        # symmetrical tests first
        expect = {'PLUS-1': (1800,2200),
                  'PLUS-3': (1800,2200),
                  'MINUS-1': (800, 1200),
                  'MINUS-3': (800, 1200)}
        for gene in genes:
            region = gene.getGene3PrimeWindowCoords(200,200)
            self.assertEqual(region, expect[gene.ensembl_id])

        # now 200 upstream, 100 downstream
        expect = {'PLUS-1': (1800,2100),
                  'PLUS-3': (1800,2100),
                  'MINUS-1': (900, 1200),
                  'MINUS-3': (900, 1200)}
        for gene in genes:
            region = gene.getGene3PrimeWindowCoords(200,100)
            self.assertEqual(region, expect[gene.ensembl_id])


class TestExpression(TestDbBase):
    reffiles = [('file-1.txt', today),
                ('file-2.txt', today)]
    
    samples = [('sample 1', 'fake sample 1'),
               ('sample 2', 'fake sample 2')]
    
    proccessed = False
    
    def setUp(self):
        """docstring for add_files_samples"""
        super(TestExpression, self).setUp()
        
        if not self.proccessed:
            add_all_gene_exons(self.session, TestGene.genes)
        
        data = [ReferenceFile(*r) for r in self.reffiles]
        data += [Sample(*s) for s in self.samples]
        self.session.add_all(data)
        self.session.commit()
        self.proccessed = True

    def test_get_sample_entries(self):
        """ make sure that samples are returned correctly """
        samples = get_sample_entries(self.session)
        names = [s.name for s in samples]
        self.assertEqual(names, ['sample 1', 'sample 2'])

    def test_get_sample_counts(self):
        """ make sure that samples are counted correctly """
        num_samples = get_sample_counts(self.session)
        self.assertEqual(num_samples, 2)
    
    def test_unique_constraint_expression(self):
        """expression records unique by probeset and reference file"""
        gene = self.session.query(Gene).filter_by(ensembl_id='PLUS-1').one()
        sample = self.session.query(Sample).filter_by(name='sample 1').one()
        reffile = self.session.query(ReferenceFile).filter_by(name='file-1.txt').all()
        reffile = reffile[0]
        data = []
        # adding multiple copies with same reffile and transcript
        for probesets, scores in [((1024, 1026), (12.3,)), ((1024, 1026), (12.3,))]:
            expressed = Expression(probesets, scores)
            expressed.reffile_id = reffile.reffile_id
            expressed.sample = sample
            expressed.gene = gene
            data.append(expressed)
        
        self.session.add_all(data)
        self.assertRaises(IntegrityError, self.session.commit)
    
    def test_unique_constraint_expressiondiff(self):
        """expression diff records unique by gene and reffile"""
        gene = self.session.query(Gene).filter_by(ensembl_id='PLUS-1').one()
        
        samples = self.session.query(Sample).all()
        reffiles = self.session.query(ReferenceFile).all()
        
        values = ((12345,), (13,), 0.1, 0)
        ediffs = []
        # now add 2 copies of one expressiondiff
        for i in range(2):
            ed = ExpressionDiff(*values)
            ed.sample_a = samples[0]
            ed.sample_b = samples[1]
            ed.reference_file = reffiles[0]
            ed.gene = gene
            ediffs.append(ed)
        
        self.session.add_all(ediffs)
        self.assertRaises(IntegrityError, self.session.commit)
    

class TestTargetGene(TestDbBase):
    reffiles = [('file-1.txt', today),
                ('file-2.txt', today)]
    
    samples = [('sample 1', 'fake sample 1'),
               ('sample 2', 'fake sample 2')]
    
    proccessed=False
    
    def setUp(self):
        """docstring for add_files_samples"""
        super(TestTargetGene, self).setUp()
        
        if not self.proccessed:
            add_all_gene_exons(self.session, TestGene.genes)
        
        data = [ReferenceFile(*r) for r in self.reffiles]
        data += [Sample(*s) for s in self.samples]
        self.session.add_all(data)
        self.session.commit()
        self.proccessed = True
    
    def test_unique_constraint_target(self):
        """study target genes can only map to single genes"""
        gene = self.session.query(Gene).filter_by(ensembl_id='PLUS-1').one()
        sample = self.session.query(Sample).filter_by(name='sample 1').one()
        reffile = self.session.query(ReferenceFile).filter_by(name='file-1.txt').one()
        data = []
        # adding multiple copies with same reffile and gene
        for i in range(2):
            e = TargetGene()
            e.reference_file = reffile
            e.gene = gene
            e.sample = sample
            data.append(e)
        
        self.session.add_all(data)
        self.assertRaises(IntegrityError, self.session.commit)

class TestQueryFunctions(TestDbBase):
    """test the db querying functions"""
    reffiles = [('file-1.txt', today),
                ('file-2.txt', today)]
    
    samples = [('sample 1', 'fake sample 1'),
               ('sample 2', 'fake sample 2')]
    
    proccessed = False
    
    def populate_db(self, **kwargs):
        singleton = kwargs.get('singleton', False)
        data = [ReferenceFile(*r) for r in self.reffiles]
        data += [Sample(*s) for s in self.samples]
        self.session.add_all(data)
        genes = self.session.query(Gene).all()
        if singleton:
            samples = self.session.query(Sample).filter_by(name='sample 1').all()
            reffiles = self.session.query(ReferenceFile).filter_by(name='file-1.txt').all()
        else:
            samples = self.session.query(Sample).all()
            reffiles = self.session.query(ReferenceFile).all()
        
        # adding multiple copies with same reffile and transcript
        for sample, reffile in zip(samples, reffiles):
            for i, gene in enumerate(genes):
                probeset, score = (1024, 21.0+i)
                expressed = Expression((probeset+i,), (score,))
                expressed.reffile_id = reffile.reffile_id
                expressed.sample = sample
                expressed.gene = gene
                self.session.add(expressed)
                
        # add a file with nothing related to it
        reffile = ReferenceFile('file-no-data.txt', today)
        self.session.add(reffile)
        self.session.commit()
    
    def setUp(self, force=False, singleton=False):
        """docstring for add_files_samples"""
        super(TestQueryFunctions, self).setUp()
        
        if not self.proccessed or force:
            add_all_gene_exons(self.session, TestGene.genes)
        
        self.populate_db(singleton=singleton)
        self.proccessed = True

    def _build_target_gene(self, reffile, gene_name, target_sample_name):
        """ helper method to simplify code for creating new TargetGene instances
        """
        tg = TargetGene()
        tg.reference_file = self.session.query(ReferenceFile).\
                filter_by(name=reffile).one()
        tg.gene = self.session.query(Gene).\
                filter_by(ensembl_id=gene_name).one()
        tg.sample = self.session.query(Sample).\
                filter_by(name=target_sample_name).one()
        return tg

    def populate_target_data(self):
        """ populates db for inclusion/exclusion tests in get_ranked_expression
        """
        # Create an extra gene which will be used by TargetGene but not in Sample
        extra_gene_dict = dict(gene=dict(ensembl_id='TARGET-1',
            symbol='agene', biotype='protein_coding', status='fake',
            description='a fake gene',
            chrom='5', start=3000, end=5000, strand=1),
            exons=[dict(ensembl_id='exon-3', rank=3, start=3050, end=4400)]
        )
        data = [Gene(**extra_gene_dict['gene'])]
        self.session.add_all(data)
        self.session.commit()

        # Create Target reference file
        reffile1 = ReferenceFile('target1.txt', today)
        reffile2 = ReferenceFile('target2.txt', today)
        reffile3 = ReferenceFile('target3.txt', today)
        data=[reffile1, reffile2, reffile3]
        self.session.add_all(data)
        self.session.commit()

        ### Create Target samples
        # target1 = Test 1 overlap, 4 sample and 2 target genes (1 match)
        target_sample1 = Sample('target 1', 'fake target 1')
        # target2 = Test 4 overlap, 4 sample and 4 target genes (4 matches)
        target_sample2 = Sample('target 2', 'fake target 2')
        # target3 = Test 0 overlap, 4 sample and 1 target gene (0 matches)
        target_sample3 = Sample('target 3', 'fake target 3')
        targets = [target_sample1, target_sample2, target_sample3]
        self.session.add_all(targets)
        self.session.commit()

        # Create one matching and one non-matching TargetGenes for target 1
        t1 = self._build_target_gene('target1.txt', 'PLUS-1', 'target 1')
        t2 = self._build_target_gene('target1.txt', 'TARGET-1', 'target 1')
        data = [t1, t2]
        self.session.add_all(data)
        self.session.commit()

        # Create four matching TargetGenes for target 2
        t1 = self._build_target_gene('target2.txt', 'PLUS-1', 'target 2')
        t2 = self._build_target_gene('target2.txt', 'PLUS-3', 'target 2')
        t3 = self._build_target_gene('target2.txt', 'MINUS-1', 'target 2')
        t4 = self._build_target_gene('target2.txt', 'MINUS-3', 'target 2')
        data = [t1, t2, t3, t4]
        self.session.add_all(data)
        self.session.commit()

        # Create 1 non-matching TargetGene for target 3
        t1 = self._build_target_gene('target3.txt', 'TARGET-1', 'target 3')
        self.session.add_all([t1])
        self.session.commit()
    
    def test_get_expression_counts(self):
        """ correctly return number of expressed genes for a sample """
        self.assertEqual(get_expression_counts(self.session), 8)
        # return correct number with/without filename
        self.assertEqual(get_expression_counts(self.session, 'sample 1'), 4)
        self.assertEqual(get_expression_counts(self.session, 'sample 1',
                data_path='file-1.txt'), 4)
        # return correct number if no records, no file
        self.assertEqual(get_expression_counts(self.session,
                'sample 1', data_path='file-no-data.txt'), 0)
        # return correct number if no records, wrong biotype
        self.assertEqual(get_expression_counts(self.session,
                'sample 1', biotype='miRNA'), 0)

    def test_get_expression_entries(self):
        """ correctly return the expressed genes for a sample """
        # this is the order they were created amd added to the table
        expect = ['sample 1', 'sample 1', 'sample 1', 'sample 1',
                  'sample 2', 'sample 2', 'sample 2', 'sample 2']
        expr = get_expression_entries(self.session)
        names = [str(e.sample.name) for e in expr]
        self.assertEqual(names, expect)
        # return correct entries with sample name
        expect = ['sample 1', 'sample 1', 'sample 1', 'sample 1']
        expr = get_expression_entries(self.session, 'sample 1')
        names = [str(e.sample.name) for e in expr]
        self.assertEqual(names, expect)
        # with file name
        expr = get_expression_entries(self.session, 'sample 1',
                data_path='file-1.txt')
        names = [str(e.sample.name) for e in expr]
        self.assertEqual(names, expect)

        # return correct entries if no records, no file
        self.assertEqual(get_expression_entries(self.session,
            'sample 1', data_path='file-no-data.txt'), [])
        # return correct entries if no records, wrong biotype
        self.assertEqual(get_expression_entries(self.session,
            'sample 1', biotype='miRNA'), [])

    def test_get_targetgene_counts(self):
        """ the number of TargetGene entries should be correct """
        self.populate_target_data()
        self.assertEqual(get_targetgene_counts(self.session), 7)
        self.assertEqual(get_targetgene_counts(self.session, 'target 1'), 2)
        self.assertEqual(get_targetgene_counts(self.session, 'target 2'), 4)
        self.assertEqual(get_targetgene_counts(self.session, 'target 3'), 1)

    def test_get_targetgene_entries(self):
        """ TargetGene entries should be correctly returned """
        self.populate_target_data()
        # this is the order they were created and added to the table
        expect1 = ['PLUS-1', 'TARGET-1', 'PLUS-1', 'PLUS-3',
                   'MINUS-1', 'MINUS-3', 'TARGET-1']

        tgs = get_targetgene_entries(self.session)
        tg_list = [str(tg.gene.ensembl_id) for tg in tgs]
        self.assertEqual(tg_list, expect1)

        expect2 = ['PLUS-1', 'TARGET-1']
        tgs = get_targetgene_entries(self.session, 'target 1')
        tg_list = [str(tg.gene.ensembl_id) for tg in tgs]
        self.assertEqual(tg_list, expect2)

    def test_get_expressed_genes_from_chrom_order(self):
        """should return the correct number of expressed genes from a chrom"""
        # Add chroms to compare against
        add_chroms(self.session, TestChrom.species, TestChrom.chromlist)

        ranked = get_genes_by_ranked_expr(self.session,
            'sample 1', chrom='2')
        for i in range(1, len(ranked)):
            self.assertTrue(ranked[i-1].Rank < ranked[i].Rank)
        
        for gene in ranked:
            self.assertTrue(gene.chrom == '2')
    
    def test_get_ranks_scores(self):
        """return correct gene mean ranks and mean scores"""
        self.setUp(force=True, singleton=True)
        genes = get_genes_by_ranked_expr(self.session,
            'sample 1')
        expected_ranks = {'PLUS-1':4, 'PLUS-3':3, 'MINUS-1':2, 'MINUS-3':1}
        expected_scores = {'PLUS-1':21, 'PLUS-3':22, 'MINUS-1':23, 'MINUS-3':24}
        for gene in genes:
            self.assertEqual(gene.Rank, expected_ranks[gene.ensembl_id])
            self.assertEqual(gene.MeanScore, expected_scores[gene.ensembl_id])
    
    def test_get_ranked_genes_order(self):
        """return correct gene order"""
        self.setUp(force=True, singleton=True)
        ranked = get_genes_by_ranked_expr(self.session, 'sample 1')
        for i in range(1, len(ranked)):
            self.assertTrue(ranked[i-1].Rank < ranked[i].Rank)
    
    def test_get_gene_entries(self):
        """ return correct genes for a release """

        genes = get_gene_entries(self.session) # returns all genes
        expected = ['PLUS-1', 'PLUS-3', 'MINUS-1', 'MINUS-3']
        gene_ids = [gene.ensembl_id for gene in genes]
        self.assertEqual(gene_ids, expected)

        genes = get_gene_entries(self.session, chrom='2') # returns chrom2 genes
        expected = ['PLUS-3', 'MINUS-1']
        gene_ids = [gene.ensembl_id for gene in genes]
        self.assertEqual(gene_ids, expected)

        genes = get_gene_entries(self.session, biotype='miRNA') # returns none
        self.assertEqual(len(genes), 0)

    def test_get_gene_counts(self):
        """ return correct counts for genes """
        self.assertEqual(get_gene_counts(self.session), 4)

        self.assertEqual(get_gene_counts(self.session, chrom='2'), 2)

        # returns none
        num_genes = get_gene_counts(self.session, biotype='miRNA')
        self.assertEqual(num_genes, 0)

    def test_query_expressed_genes_with_inclusive_target_genes(self):
        """ return only those genes which overlap with target """
        # Need to test when we have 100% overlap, 0% overlap and something in between

        self.populate_target_data()

        # Test 1 overlap, 4 sample and 2 target genes
        included_genes = get_genes_by_ranked_expr(self.session, 'sample 1',
                include_targets='target 1')
        self.assertTrue(len(included_genes) == 1)
        self.assertTrue(included_genes[0].ensembl_id == 'PLUS-1')

        # Test 4 overlap, 4 sample and 4 target genes
        included_genes = get_genes_by_ranked_expr(self.session, 'sample 1',
                include_targets='target 2')
        self.assertTrue(len(included_genes) == 4)
        self.assertTrue(included_genes[0].ensembl_id == 'MINUS-3')
        self.assertTrue(included_genes[1].ensembl_id == 'MINUS-1')
        self.assertTrue(included_genes[2].ensembl_id == 'PLUS-3')
        self.assertTrue(included_genes[3].ensembl_id == 'PLUS-1')

        # Test 0 overlap, 4 sample and 1 non-matching target gene
        remaining_genes = get_genes_by_ranked_expr(self.session, 'sample 1',
                include_targets='target 3')
        self.assertTrue(len(remaining_genes) == 0)

    def test_query_expressed_genes_with_exclusive_target_genes(self):
        """ return only those genes which DON'T overlap with target """
        # Need to test when we have 100% overlap, 0% overlap and something in between

        self.populate_target_data()

        # Test 1 overlap, 4 sample and 2 target genes
        non_overlapping_genes = get_genes_by_ranked_expr(self.session,
                'sample 1', exclude_targets='target 1')
        self.assertTrue(len(non_overlapping_genes) == 3)
        self.assertTrue(non_overlapping_genes[0].ensembl_id == 'MINUS-3')
        self.assertTrue(non_overlapping_genes[1].ensembl_id == 'MINUS-1')
        self.assertTrue(non_overlapping_genes[2].ensembl_id == 'PLUS-3')

        # Test 4 overlap, 4 sample and 4 target genes
        non_overlapping_genes = get_genes_by_ranked_expr(self.session, 'sample 1',
            exclude_targets='target 2')
        self.assertTrue(len(non_overlapping_genes) == 0)

        # Test 0 overlap, 4 sample and 1 non-matching target gene
        non_overlapping_genes = get_genes_by_ranked_expr(self.session, 'sample 1',
            exclude_targets='target 3')
        self.assertTrue(len(non_overlapping_genes) == 4)
        self.assertTrue(non_overlapping_genes[0].ensembl_id == 'MINUS-3')
        self.assertTrue(non_overlapping_genes[1].ensembl_id == 'MINUS-1')
        self.assertTrue(non_overlapping_genes[2].ensembl_id == 'PLUS-3')
        self.assertTrue(non_overlapping_genes[3].ensembl_id == 'PLUS-1')


class TestQueryFunctionsExpDiff(TestDbBase):
    """test the db querying functions"""
    reffiles = [('file-1.txt', today),
                ('file-2.txt', today)]
    
    samples = [('sample 1', 'fake sample 1'),
               ('sample 2', 'fake sample 2')]
    
    dpath = 'data/microarray_diff.txt'
    sample = ('sample1', 'blah')
    reffile_path1 = 'sample1.txt'
    reffile_path2 = 'sample2.txt'
    
    def populate_db(self, **kwargs):
        # setting up some starting values
        reffile1 = ReferenceFile(self.reffile_path1, today)
        reffile2 = ReferenceFile(self.reffile_path2, today)
        name, desc = self.sample
        success = add_sample(self.session, name, desc)
        self.assertTrue(success)
        self.session.add_all([reffile1, reffile2])
        self.session.commit()
        table = gene_expr_diff_to_table(self.dpath, stable_id_label='gene',
                probeset_label='probeset', exp_label='exp', sig_label='sig',
                pval_label='p_val')
        add_expression_diff_study(self.session, 'sample1', self.dpath, table,
                self.reffile_path1, self.reffile_path2, show_progress=False)
    
    def setUp(self):
        """docstring for add_files_samples"""
        super(TestQueryFunctionsExpDiff, self).setUp()
        add_all_gene_exons(self.session, TestGene.genes)
        self.populate_db()

    def _build_target_gene(self, reffile, gene_name, target_sample_name):
        """ helper method to simplify code for creating new TargetGene instances
        """
        tg = TargetGene()
        tg.reference_file = self.session.query(ReferenceFile).\
        filter_by(name=reffile).one()
        tg.gene = self.session.query(Gene).\
        filter_by(ensembl_id=gene_name).one()
        tg.sample = self.session.query(Sample).\
        filter_by(name=target_sample_name).one()
        return tg

    def populate_target_data(self):
        """ populates db for inclusion/exclusion tests in get_ranked_expression
        """
        # Create an extra gene which will be used by TargetGene but not in Sample
        extra_gene_dict = dict(gene=dict(ensembl_id='TARGET-1',
            symbol='agene', biotype='protein_coding', status='fake',
            description='a fake gene',
            chrom='5', start=3000, end=5000, strand=1),
            exons=[dict(ensembl_id='exon-3', rank=3, start=3050, end=4400)]
        )
        data = [Gene(**extra_gene_dict['gene'])]
        self.session.add_all(data)
        self.session.commit()

        # Create Target reference file
        reffile1 = ReferenceFile('target1.txt', today)
        reffile2 = ReferenceFile('target2.txt', today)
        reffile3 = ReferenceFile('target3.txt', today)
        data=[reffile1, reffile2, reffile3]
        self.session.add_all(data)
        self.session.commit()

        ### Create Target samples
        # target1 = Test 1 overlap, 4 sample and 2 target genes (1 match)
        target_sample1 = Sample('target 1', 'fake target 1')
        # target2 = Test 4 overlap, 4 sample and 4 target genes (4 matches)
        target_sample2 = Sample('target 2', 'fake target 2')
        # target3 = Test 0 overlap, 4 sample and 1 target gene (0 matches)
        target_sample3 = Sample('target 3', 'fake target 3')
        targets = [target_sample1, target_sample2, target_sample3]
        self.session.add_all(targets)
        self.session.commit()

        # Create one matching and one non-matching TargetGenes for target 1
        t1 = self._build_target_gene('target1.txt', 'PLUS-1', 'target 1')
        t2 = self._build_target_gene('target1.txt', 'TARGET-1', 'target 1')
        data = [t1, t2]
        self.session.add_all(data)
        self.session.commit()

        # Create four matching TargetGenes for target 2
        t1 = self._build_target_gene('target2.txt', 'PLUS-1', 'target 2')
        t2 = self._build_target_gene('target2.txt', 'PLUS-3', 'target 2')
        t3 = self._build_target_gene('target2.txt', 'MINUS-1', 'target 2')
        t4 = self._build_target_gene('target2.txt', 'MINUS-3', 'target 2')
        data = [t1, t2, t3, t4]
        self.session.add_all(data)
        self.session.commit()

        # Create 1 non-matching TargetGene for target 3
        t1 = self._build_target_gene('target3.txt', 'TARGET-1', 'target 3')
    
    def test_add_expression_diff_data(self):
        """correctly add expression difference data"""
        # add the expression diff data
        # do we get it back?
        expr_diffs = get_diff_entries(self.session, self.sample[0])
        expect = dict([('PLUS-1', ['10600707']),
                  ('PLUS-3', ['10408081']),
                  ('MINUS-1', ['10494402']),
                  ('MINUS-3', ['10408083'])])

        self.assertTrue(len(expr_diffs) > 0)
        for diff in expr_diffs:
            expect_probeset = expect[diff.gene.ensembl_id]
            self.assertEqual(diff.probesets, expect_probeset)
    
    def test_query_exp_diff(self):
        """return correct records from query when filtered"""
        name_start = {-1: 'MINUS', 1: 'PLUS'}
        for multitest_signif_val in [-1, 1]:
            expr_diffs = get_diff_entries(self.session, self.sample[0],
                multitest_signif_val=multitest_signif_val)
            
            self.assertEqual(len(expr_diffs), 2)
            # should only get records with the correct test significance
            # direction
            for record in expr_diffs:
                self.assertEqual(record.multitest_signif,
                    multitest_signif_val)
                # gene names designed to match the test significance
                self.assertTrue(record.gene.ensembl_id.startswith(
                    name_start[multitest_signif_val]))

    def test_get_expressed_diff_genes_from_chrom_order(self):
        """should return the correct number of expressed genes from a chrom"""
        # Add chroms to compare against
        add_chroms(self.session, TestChrom.species, TestChrom.chromlist)

        ranked = get_genes_by_ranked_diff(self.session,
                'sample1', multitest_signif_val=1, chrom='2')
        for i in range(1, len(ranked)):
            self.assertTrue(ranked[i-1].Rank < ranked[i].Rank)

        for gene in ranked:
            self.assertTrue(gene.chrom == '2')
    
    def test_get_expressed_diff_genes(self):
        """return genes ranked by foldchange"""
        genes = get_genes_by_ranked_diff(self.session, self.sample[0])
        self.assertTrue(len(genes) == 4)
        for i in range(3):
            self.assertTrue(genes[i].Rank < genes[i+1].Rank)
        
        # sample up genes
        genes = get_genes_by_ranked_diff(self.session, self.sample[0],
            multitest_signif_val=1)
        self.assertTrue(len(genes) == 2)
        expect_order = ['PLUS-1', 'PLUS-3']
        for i in range(2):
            self.assertEqual(genes[i].ensembl_id, expect_order[i])
        
        # sample down genes
        genes = get_genes_by_ranked_diff(self.session, self.sample[0],
            multitest_signif_val=-1)
        self.assertTrue(len(genes) == 2)
        expect_order = ['MINUS-1', 'MINUS-3']
        for i in range(2):
            self.assertEqual(genes[i].ensembl_id, expect_order[i])

    def test_query_expressed_diff_genes_with_inclusive_target_genes(self):
        """ return only those genes which overlap with target """
        # Need to test when we have 100% overlap, 0% overlap and something in between

        self.populate_target_data()

        # Test 1 overlap, 4 sample and 2 target genes
        remaining_genes = get_genes_by_ranked_diff(self.session, 'sample1',
            include_targets='target 1')
        self.assertTrue(len(remaining_genes) == 1)
        self.assertTrue(remaining_genes[0].ensembl_id == 'PLUS-1')

        # Test 4 overlap, 4 sample and 4 target genes
        remaining_genes = get_genes_by_ranked_diff(self.session, 'sample1',
            include_targets='target 2')
        self.assertTrue(len(remaining_genes) == 4)
        self.assertTrue(remaining_genes[0].ensembl_id == 'PLUS-1')
        self.assertTrue(remaining_genes[1].ensembl_id == 'PLUS-3')
        self.assertTrue(remaining_genes[2].ensembl_id == 'MINUS-1')
        self.assertTrue(remaining_genes[3].ensembl_id == 'MINUS-3')

        # Test 0 overlap, 4 sample and 1 non-matching target gene
        remaining_genes = get_genes_by_ranked_diff(self.session, 'sample1',
            include_targets='target 3')
        self.assertTrue(len(remaining_genes) == 0)

    def test_query_expressed_diff_genes_with_exclusive_target_genes(self):
        """ return only those genes which DON'T overlap with target """
        # Need to test when we have 100% overlap, 0% overlap and something in between

        self.populate_target_data()

        # Test 1 overlap, 4 sample and 2 target genes
        remaining_genes = get_genes_by_ranked_diff(self.session, 'sample1',
            exclude_targets='target 1')
        self.assertTrue(len(remaining_genes) == 3)
        self.assertTrue(remaining_genes[0].ensembl_id == 'PLUS-3')
        self.assertTrue(remaining_genes[1].ensembl_id == 'MINUS-1')
        self.assertTrue(remaining_genes[2].ensembl_id == 'MINUS-3')

        # Test 4 overlap, 4 sample and 4 target genes
        remaining_genes = get_genes_by_ranked_diff(self.session, 'sample1',
            exclude_targets='target 2')
        self.assertTrue(len(remaining_genes) == 0)

        # Test 0 overlap, 4 sample and 1 non-matching target gene
        remaining_genes = get_genes_by_ranked_diff(self.session, 'sample1',
            exclude_targets='target 3')
        self.assertTrue(len(remaining_genes) == 4)
        self.assertTrue(remaining_genes[0].ensembl_id == 'PLUS-1')
        self.assertTrue(remaining_genes[1].ensembl_id == 'PLUS-3')
        self.assertTrue(remaining_genes[2].ensembl_id == 'MINUS-1')
        self.assertTrue(remaining_genes[3].ensembl_id == 'MINUS-3')

    def test_get_diff_entries(self):
        """ make sure diff table entries match expected """
        entries = get_diff_entries(self.session)
        self.assertEqual(entries[0].expression_diff_id, 1)
        self.assertEqual(entries[1].expression_diff_id, 2)
        self.assertEqual(entries[2].expression_diff_id, 3)
        self.assertEqual(entries[3].expression_diff_id, 4)

    def test_get_diff_counts(self):
        """ make sure diff table counts match expected """
        self.assertEqual(get_diff_counts(self.session), 4)

    def test_drop_sample_records(self):
        """ test the removal of all DB records based on sample name """

        self.populate_target_data()

        expected = get_sample_counts(self.session)
        # Try with built-in test mode, which rolls back the deletes
        self.assertFalse(drop_sample_records(self.session,
                'sample1', test=True))
        # Now actually delete
        self.assertTrue(drop_sample_records(self.session,
                'sample1', test=False))

        observed = get_sample_counts(self.session)
        self.assertEqual(observed, expected-1)

if __name__ == '__main__':
    main()
