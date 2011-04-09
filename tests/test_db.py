"""test design of the express.db module"""
import sys
sys.path.append('..')

import warnings
warnings.filterwarnings('ignore', 'Not using MPI as mpi4py not found')

from cogent.util.unit_test import TestCase, main

import datetime
from sqlalchemy import create_engine, and_, or_
from sqlalchemy.exc import IntegrityError

from chippy.express.db_schema import Gene, Exon, \
            ExternalGene, Expression, ExpressionDiff, ReferenceFile, Sample, \
            make_session
from chippy.express.db_query import get_total_gene_counts, \
        get_ranked_expression, get_ranked_genes_per_chrom

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
        self.session = make_session("sqlite:///:memory:")

ensembl_release = '58'

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
        reffiles = self.session.query(ReferenceFile).all()
        self.assertEqual(reffiles[0].sample.name, 'A')


class TestGene(TestDbBase):
    """test gene properties"""
    plus_coords_one_exons = dict(gene=dict(ensembl_id='PLUS-1',
        ensembl_release=ensembl_release,
        symbol='agene', biotype='protein_coding', status='fake',
        description='a fake gene',
        coord_name='1', start=1000, end=2000, strand=1),
        exons=[dict(ensembl_id='exon-1', rank=1, start=1050, end=1950,
                    ensembl_release=ensembl_release)]
        )
    
    plus_coords_many_exons = dict(gene=dict(ensembl_id='PLUS-3',
        ensembl_release='58',
        symbol='agene', biotype='protein_coding', status='fake',
        description='a fake gene',
        coord_name='2', start=1000, end=2000, strand=1), 
        exons=[dict(ensembl_id='exon-1', rank=1, start=1050, end=1400,
                    ensembl_release=ensembl_release),
               dict(ensembl_id='exon-2', rank=2, start=1600, end=1700,
                    ensembl_release=ensembl_release),
               dict(ensembl_id='exon-3', rank=3, start=1800, end=1900,
                    ensembl_release=ensembl_release)]
        )
    
    # 
    minus_coords_one_exons = dict(gene=dict(ensembl_id='MINUS-1',
        ensembl_release=ensembl_release,
        symbol='agene', biotype='protein_coding', status='fake',
        description='a fake gene',
        coord_name='2', start=1000, end=2000, strand=-1),
        exons=[dict(ensembl_id='exon-1', rank=1, start=1050, end=1950,
                    ensembl_release=ensembl_release)]
        )
    
    minus_coords_many_exons = dict(gene=dict(ensembl_id='MINUS-3',
        ensembl_release=ensembl_release,
        symbol='agene', biotype='protein_coding', status='fake',
        description='a fake gene',
        coord_name='3', start=1000, end=2000, strand=-1), 
        exons=[dict(ensembl_id='exon-3', rank=3, start=1050, end=1400,
                    ensembl_release=ensembl_release),
               dict(ensembl_id='exon-2', rank=2, start=1600, end=1700,
                    ensembl_release=ensembl_release),
               dict(ensembl_id='exon-1', rank=1, start=1800, end=1900,
                    ensembl_release=ensembl_release)]
        )
    
    genes = [plus_coords_one_exons, plus_coords_many_exons,
            minus_coords_one_exons, minus_coords_many_exons]
    
    def test_add_genes(self):
        """excercise adding a gene"""
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
    
    def test_gene_tss(self):
        """return correct TSS coordinate"""
        genes = self.session.query(Gene).all()
        expect = {'PLUS-1': 1000,
            'PLUS-3': 1000,
            'MINUS-1': 2000,
            'MINUS-3': 2000}
        for gene in genes:
            self.assertEqual(gene.Tss, expect[gene.ensembl_id])
    
    def test_tss_centred_coords(self):
        """return coordinates centred on the TSS"""
        add_all_gene_exons(self.session, self.genes)
        genes = self.session.query(Gene).all()
        expect = {'PLUS-1': [500, 1500],
            'PLUS-3': [500, 1500],
            'MINUS-1': [1500, 2500],
            'MINUS-3': [1500, 2500]}
        for gene in genes:
            self.assertEqual(gene.getTssCentredCoords(500),
                            expect[gene.ensembl_id])
    
    def test_gene_upstream(self):
        """return correct coordinates ending at TSS"""
        add_all_gene_exons(self.session, self.genes)
        genes = self.session.query(Gene).all()
        expect = {'PLUS-1': (500,1000),
            'PLUS-3': (500, 1000),
            'MINUS-1': (2000, 2500),
            'MINUS-3': (2000, 2500)}
        for gene in genes:
            self.assertEqual(gene.getUpstreamCoords(500),
                            expect[gene.ensembl_id])
    

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
    
    def est_unique_constraint_expressiondiff(self):
        """expression diff records unique by id of expression in samples A & B"""
        # TODO turned this test off since there is a refactor taking place
        # and the schema is not yet finalised for this case
        gene = self.session.query(Gene).filter_by(ensembl_id='PLUS-1').one()
        
        samples = self.session.query(Sample).all()
        reffiles = self.session.query(ReferenceFile).all()
        
        values = (12345, 13, 0.1, 0)
        ediffs = []
        # now add 2 copies of one expressiondiff
        for i in range(2):
            ed = ExpressionDiff(*values)
            ed.sample_a = samples[0]
            ed.sample_b = samples[1]
            ed.reference_file = reffiles[0]
            ed.transcript = canonical_transcript
            ediffs.append(ed)
        
        self.session.add_all(ediffs)
        self.assertRaises(IntegrityError, self.session.commit)
    

class TestExternalGene(TestDbBase):
    reffiles = [('file-1.txt', today),
                ('file-2.txt', today)]
    
    samples = [('sample 1', 'fake sample 1'),
               ('sample 2', 'fake sample 2')]
    
    proccessed=False
    
    def setUp(self):
        """docstring for add_files_samples"""
        super(TestExternalGene, self).setUp()
        
        if not self.proccessed:
            add_all_gene_exons(self.session, TestGene.genes)
        
        data = [ReferenceFile(*r) for r in self.reffiles]
        data += [Sample(*s) for s in self.samples]
        self.session.add_all(data)
        self.session.commit()
        self.proccessed = True
    
    def test_unique_constraint_external(self):
        """study external genes can only map to single genes"""
        gene = self.session.query(Gene).filter_by(ensembl_id='PLUS-1').one()
        sample = self.session.query(Sample).filter_by(name='sample 1').one()
        reffile = self.session.query(ReferenceFile).filter_by(name='file-1.txt').one()
        data = []
        # adding multiple copies with same reffile and gene
        for i in range(2):
            e = ExternalGene()
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
                
        # add a file with nothign related to it
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
    
    def test_counting_genes(self):
        """correctly return number of genes for a sample"""
        # return correct number with/without filename
        self.assertEqual(get_total_gene_counts(self.session, ensembl_release,
            'sample 1'), 4)
        self.assertEqual(get_total_gene_counts(self.session, ensembl_release,
            'sample 1', data_path='file-1.txt'), 4)
        # return correct number if no records, no file
        self.assertEqual(get_total_gene_counts(self.session, ensembl_release,
            'sample 1', data_path='file-no-data.txt'), 0)
        # return correct number if no records, wrong biotype
        self.assertEqual(get_total_gene_counts(self.session, ensembl_release,
            'sample 1', biotype='miRNA'), 0)
    
    def test_get_expressed_genes_from_chrom(self):
        """should return the correct number of expressed genes from a chrom"""
        ranked = get_ranked_genes_per_chrom(self.session, ensembl_release,
            'sample 1', '2')
        for i in range(1, len(ranked)):
            self.assertTrue(ranked[i-1].Rank < ranked[i].Rank)
        
        for gene in ranked:
            self.assertTrue(gene.coord_name == '2')
    
    def test_get_ranks_scores(self):
        """return correct gene mean ranks and mean scores"""
        self.setUp(force=True, singleton=True)
        genes = get_ranked_expression(self.session, ensembl_release,
            'sample 1')
        expected_ranks = {'PLUS-1':4, 'PLUS-3':3, 'MINUS-1':2, 'MINUS-3':1}
        expected_scores = {'PLUS-1':21, 'PLUS-3':22, 'MINUS-1':23, 'MINUS-3':24}
        for gene in genes:
            self.assertEqual(gene.Rank,
                        expected_ranks[gene.ensembl_id])
            self.assertEqual(gene.MeanScore,
                        expected_scores[gene.ensembl_id])
    
    def test_get_ranked_genes(self):
        """return correct gene order"""
        self.setUp(force=True, singleton=True)
        ranked = get_ranked_expression(self.session, ensembl_release,
            'sample 1')
        for i in range(1, len(ranked)):
            self.assertTrue(ranked[i-1].Rank < ranked[i].Rank)
        
    

if __name__ == '__main__':
    main()
