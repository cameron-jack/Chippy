#!/usr/bin/env python
from __future__ import division

import os, sys, glob
sys.path.extend(['../../src'])

from optparse import make_option

from cogent import LoadTable
from cogent.util.progress_display import display_wrap
from cogent.util.misc import parse_command_line_parameters

from chippy.express import db_query, db_schema
from chippy.core.read_count import get_combined_counts
from chippy.ref.util import chroms
from chippy.util.util import grouped_by_chrom
from chippy.util.run_record import RunRecord
from chippy.util.definition import LOG_DEBUG, LOG_INFO, LOG_WARNING, \
    LOG_ERROR, LOG_CRITICAL

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2011"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "0.9.dev"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Development"

if 'CHIPPY_DB' in os.environ:
    db_path = os.environ['CHIPPY_DB']
else:
    raise RuntimeError('You need to set an environment variable CHIPPY_DB '\
                       'that indicates where to find the database')

session = db_query.make_session('sqlite:///%s' % db_path)
samples = db_query.get_samples(session)
if not samples:
    samples = [None]

def _get_exon_counts(gene, counts):
    rows = []
    for rank, (start, end) in enumerate(gene.ExonCoordsByRank):
        size = (end-start)
        tag_counts = counts[start:end].sum()
        rows.append(['exon', gene.ensembl_id, rank, tag_counts, size])
    return rows

def _get_intron_counts(gene, counts):
    rows = []
    for rank, (start, end) in enumerate(gene.IntronCoordsByRank):
        size = (end-start)
        tag_counts = counts[start:end].sum()
        rows.append(['intron', gene.ensembl_id, rank, tag_counts, size])
    return rows

def _get_upstream_counts(gene, counts, size):
    rows = []
    start, end = gene.getUpstreamCoords(size)
    size = (end-start)
    tag_counts = counts[start:end].sum()
    rows.append(['upstream', gene.ensembl_id, 0, tag_counts, size])
    return rows

def _get_count_sum_table_per_chrom(counts, genes, upstream_size):
    """returns table of total counts for upstream, exon, intron coords
    """
    rows = []
    header = ['region_type', 'ensembl_id', 'region_rank', 'counts', 'size']
    
    for gene in genes:
        # if no intron, we discard
        if len(gene.IntronCoords) == 0:
            continue
        
        rows += _get_exon_counts(gene, counts)
        rows += _get_intron_counts(gene, counts)
        rows += _get_upstream_counts(gene, counts, upstream_size)
    
    table = LoadTable(header=header, rows=rows)
    return table

def get_sum_counts_table(session, chrom_genes, counts_dir, max_read_length, count_max_length, upstream_size, test_run=False):
    """returns summary counts table"""
    dirname = os.path.dirname(counts_dir)
    basename = os.path.basename(counts_dir)
    counts_dirs = [os.path.join(dirname, p) for p in glob.glob1(dirname,
                                                    basename)]
    if test_run: # do smaller chromosomes only
        chrom_names = ['18', '19']
    else:
        chrom_names = sorted([(len(n), n) for n in chrom_genes])
        chrom_names = [n for l,n in chrom_names]
    
    tables = []
    for chrom_name in chrom_names:
        chrom_counts = get_combined_counts(counts_dirs, chrom_name,
                    max_read_length, count_max_length)
        chrom_table = _get_count_sum_table_per_chrom(chrom_counts,
                                    chrom_genes[chrom_name], upstream_size)
        chrom_table.Title = chrom_name
        tables.append(chrom_table)
        del chrom_counts
    
    table = tables.pop(0)
    table = table.appended('coord_name', tables)
    table.Title = ''
    return table

script_info = {}

script_info['title'] = 'Exports summed counts for genome feature types'
script_info['script_description'] = "Writes a compressed tab delimited"\
    " summary table that can be used for conducting statistical analyses."\
    " NOTE: Genes without an intron are skipped."
script_info['version'] = __version__

script_info['required_options'] = [
make_option('-c', '--sample', type='choice',
           help='Choose the expression study [default: %default]',
           choices=[str(s) for s in samples]),
    make_option('-1','--IP_counts_path',
        help='path to directory containing immuno precipitated mapped read'\
             ' summary tables. Can be a glob pattern.'),
    make_option('-2','--IN_counts_path',
        help='path to directory containing input mapped read summary tables.'\
            'Can be a glob pattern.'),
    make_option('-s','--save_table_name', help='path to save to table data'),
    make_option('-S', '--sample_top', type='int', default = None,
     help='Genes ranked to this number will be sampled (default is All)'),
]

script_info['optional_options'] = [
    make_option('-u',
             '--upstream_size',
             type='int', default=500,
             help='Size of region upstream from gene'),
    make_option('-m',
             '--maximum_read_length',
             type='int', default=75,
             help='Maximum sequence read length'),
    make_option('-k', '--count_max_length',
             action='store_false',
             help="Use maximum read length instead of mapped length",
             default=True),
    make_option('-p',
             '--pseudo_count',
             type='int', default=1,
             help='Value to use in ratio if input has no counts'),
    make_option('-t',
             '--test_run',
             action='store_true',
             help='No files will be written',
             default=False),
]

script_info['authors'] = __author__

def _samples_name(sample):
    """returns name from 'name : description' string"""
    name = sample.split(':')[0].strip()
    return name

def CalcRatio(pseudo_count):
    def call((x,y)):
        y = y or pseudo_count
        return x/y
    return call

def main():
    option_parser, opts, args =\
    parse_command_line_parameters(**script_info)
    
    rr = RunRecord()
    rr.addMessage('export_feature_counts', LOG_INFO,
        'Chosen sample', opts.sample)
    
    if opts.sample is None:
        rr.display()
        raise RuntimeError('No samples available')
        return
    
    sample_name = _samples_name(opts.sample)
    rr.addMessage('export_feature_counts', LOG_INFO,
        'Chosen sample', opts.sample)
    rr.addMessage('export_feature_counts', LOG_INFO,
        'Immuno precipitated counts path', opts.IP_counts_path)
    rr.addMessage('export_feature_counts', LOG_INFO,
        'Input counts path', opts.IN_counts_path)
    rr.addMessage('export_feature_counts', LOG_INFO,
        'No. of most expressed genes sampled', opts.sample_top)
    rr.addMessage('export_feature_counts', LOG_INFO,
        'maximum_read_length', opts.maximum_read_length)
    rr.addMessage('export_feature_counts', LOG_INFO,
        'count_max_length', opts.count_max_length)
    rr.addMessage('export_feature_counts', LOG_INFO,
        'upstream_size', opts.upstream_size)
    rr.addMessage('export_feature_counts', LOG_INFO,
        'pseudo_count', opts.pseudo_count)
    
    genes = db_query.get_ranked_expression(session,
                sample_name, biotype='protein_coding', rank_by='mean',
                test_run=opts.test_run)
    genes = genes[:opts.sample_top]
    chrom_gene_groups = grouped_by_chrom(genes)
    
    ip_table = get_sum_counts_table(session, chrom_gene_groups,
            opts.IP_counts_path, opts.maximum_read_length,
            opts.count_max_length, opts.upstream_size, test_run=opts.test_run)
    # IP stands for Immuno precipitated
    ip_table.Title = 'IP'
    # renamed the IP table counts header for consistency with result of
    # joining the IP and IN tables
    ip_table = ip_table.withNewHeader(['counts'], ['IP_counts'])
    
    in_table = get_sum_counts_table(session, chrom_gene_groups,
            opts.IN_counts_path, opts.maximum_read_length,
            opts.count_max_length, opts.upstream_size, test_run=opts.test_run)
    
    # IN stands for Input
    in_table.Title = 'IN'
    
    combined = ip_table.joined(in_table,
                    columns_self=('region_type', 'ensembl_id', 'region_rank'))
    combined.Title = ''
    
    ratio = CalcRatio(opts.pseudo_count)
    combined = combined.withNewColumn('ratio', ratio,
                columns=['IP_counts', 'IN_counts'])
    
    if not opts.test_run:
        combined.writeToFile(opts.save_table_name, sep='\t')
        rr.addMessage('export_feature_counts', LOG_INFO,
            'Wrote counts to', opts.save_table_name)
    else:
        print combined[:10]
    
    rr.display()
    

if __name__ == "__main__":
    main()

