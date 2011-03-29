"""test design of the express.db module"""
import sys
sys.path.append('..')

from cogent.util.unit_test import TestCase, main

import datetime
from sqlalchemy import create_engine, and_, or_
from sqlalchemy.exc import IntegrityError

from chippy.express.db import Gene, Transcript, Exon, \
            ExternalGene, Expression, ExpressionDiff, ReferenceFile, Sample, \
            make_session

now = datetime.datetime.now()
today = datetime.date(now.year, now.month, now.day)

def add_all_gene_transcript_exons(session, genes):
    data = []
    for record in genes:
        gene_data = record['gene']
        exons_data = record['exons']
        gene = Gene(**gene_data)
        data.append(gene)
        ts = Transcript(**record['transcripts'])
        ts.gene = gene
        for e_data in exons_data:
            exon = Exon(**e_data)
            exon.gene = gene
            exon.transcript = ts
            data.append(exon)
    session.add_all(data)


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
                    ensembl_release=ensembl_release)],
        transcripts=dict(ensembl_id='PLUS-1-trans-1',
                    ensembl_release=ensembl_release)
        )
    
    plus_coords_many_exons = dict(gene=dict(ensembl_id='PLUS-3',
        ensembl_release='58',
        symbol='agene', biotype='protein_coding', status='fake',
        description='a fake gene',
        coord_name='1', start=1000, end=2000, strand=1), 
        exons=[dict(ensembl_id='exon-1', rank=1, start=1050, end=1400,
                    ensembl_release=ensembl_release),
               dict(ensembl_id='exon-2', rank=2, start=1600, end=1700,
                    ensembl_release=ensembl_release),
               dict(ensembl_id='exon-3', rank=3, start=1800, end=1900,
                    ensembl_release=ensembl_release)],
        transcripts=dict(ensembl_id='PLUS-3-trans-1',
                    ensembl_release=ensembl_release))
    
    # 
    minus_coords_one_exons = dict(gene=dict(ensembl_id='MINUS-1',
        ensembl_release=ensembl_release,
        symbol='agene', biotype='protein_coding', status='fake',
        description='a fake gene',
        coord_name='1', start=1000, end=2000, strand=-1),
        exons=[dict(ensembl_id='exon-1', rank=1, start=1050, end=1950,
                    ensembl_release=ensembl_release)],
        transcripts=dict(ensembl_id='MINUS-1-trans-1',
                    ensembl_release=ensembl_release)
        )
    
    minus_coords_many_exons = dict(gene=dict(ensembl_id='MINUS-3',
        ensembl_release=ensembl_release,
        symbol='agene', biotype='protein_coding', status='fake',
        description='a fake gene',
        coord_name='1', start=1000, end=2000, strand=-1), 
        exons=[dict(ensembl_id='exon-3', rank=1, start=1050, end=1400,
                    ensembl_release=ensembl_release),
               dict(ensembl_id='exon-2', rank=2, start=1600, end=1700,
                    ensembl_release=ensembl_release),
               dict(ensembl_id='exon-1', rank=3, start=1800, end=1900,
                    ensembl_release=ensembl_release)],
        transcripts=dict(ensembl_id='MINUS-3-trans-1',
                    ensembl_release=ensembl_release))
    
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
        ts = Transcript(**self.plus_coords_many_exons['transcripts'])
        ts.gene = gene
        data.append(gene)
        data.append(ts)
        for e_data in self.plus_coords_many_exons['exons']:
            exon = Exon(**e_data)
            exon.gene = gene
            exon.transcript = ts
            data.append(exon)
            
        for e_data in self.plus_coords_many_exons['exons']:
            exon = Exon(**e_data)
            exon.gene = gene
            exon.transcript = ts
            data.append(exon)
        
        self.session.add_all(data)
        self.assertRaises(IntegrityError, self.session.commit)
    
    def test_get_gene_exon_coords(self):
        """Gene instances correctly derive coords for their exons"""
        add_all_gene_transcript_exons(self.session, self.genes)
        genes = self.session.query(Gene).all()
        expect = {'PLUS-1': [(1050, 1950)],
            'PLUS-3': [(1050, 1400),(1600, 1700),(1800, 1900)],
            'MINUS-1': [(1050, 1950)],
            'MINUS-3': [(1050, 1400),(1600, 1700),(1800, 1900)],}
        for gene in genes:
            self.assertEqual(gene.ExonCoords, expect[gene.ensembl_id])
    
    def test_get_gene_intron_coords(self):
        """Gene instances correctly derive coords for their intron"""
        add_all_gene_transcript_exons(self.session, self.genes)
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
        add_all_gene_transcript_exons(self.session, self.genes)
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
        add_all_gene_transcript_exons(self.session, self.genes)
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
    
    
    
    proccessed=False
    
    def setUp(self):
        """docstring for add_files_samples"""
        super(TestExpression, self).setUp()
        
        if not self.proccessed:
            add_all_gene_transcript_exons(self.session, TestGene.genes)
        
        data = [ReferenceFile(*r) for r in self.reffiles]
        data += [Sample(*s) for s in self.samples]
        self.session.add_all(data)
        self.session.commit()
        self.proccessed = True
    
    def test_unique_constraint_expression(self):
        """expression records unique by gene and reference file"""
        gene = self.session.query(Gene).filter_by(ensembl_id='PLUS-1').one()
        sample = self.session.query(Sample).filter_by(name='sample 1').one()
        reffile = self.session.query(ReferenceFile).filter_by(name='file-1.txt').all()
        reffile = reffile[0]
        data = []
        # adding multiple copies with same reffile and gene
        for score, rank in [(12.3, 1), (12.3, 1)]:
            e = Expression(score, rank)
            e.reference_file = reffile
            e.gene = gene
            e.sample = sample
            data.append(e)
        
        self.session.add_all(data)
        self.assertRaises(IntegrityError, self.session.commit)
    
    def test_unique_constraint_expressiondiff(self):
        """expression diff records unique by id of expression in samples A & B"""
        gene = self.session.query(Gene).filter_by(ensembl_id='PLUS-1').one()
        samples = self.session.query(Sample).all()
        reffiles = self.session.query(ReferenceFile).all()
        
        values = (13, 0.1, 0)
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
            add_all_gene_transcript_exons(self.session, TestGene.genes)
        
        data = [ReferenceFile(*r) for r in self.reffiles]
        data += [Sample(*s) for s in self.samples]
        self.session.add_all(data)
        self.session.commit()
        self.proccessed = True
    
    def test_unique_constraint_expression(self):
        """docstring for test_a"""
        gene = self.session.query(Gene).filter_by(ensembl_id='PLUS-1').all()
        gene = gene[0]
        sample = self.session.query(Sample).filter_by(name='sample 1').all()
        sample = sample[0]
        reffile = self.session.query(ReferenceFile).filter_by(name='file-1.txt').all()
        reffile = reffile[0]
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
        


if __name__ == '__main__':
    main()
