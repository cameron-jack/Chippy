#!/usr/bin/env python
'''reports a summary of the database - needs more work'''

from __future__ import division

import sys
sys.path.extend(['..'])

from chippy.express.db_query import make_session, get_chroms, get_species,\
        get_sample_counts, get_reffile_counts, get_gene_counts, get_exon_counts,\
        get_expression_counts, get_diff_counts, get_targetgene_counts,\
        get_reffile_entries, get_all_sample_names
from chippy.util.command_args import Args
from chippy.util.run_record import RunRecord

__author__ = 'Gavin Huttley, Cameron Jack'
__copyright__ = 'Copyright 2012, Gavin Huttley, Cameron Jack, Anuj Pahwa'
__credits__ = ['Gavin Huttley', 'Cameron Jack']
__license__ = 'GPL'
__version__ = '746'
__maintainer__ = 'Cameron Jack'
__email__ = 'cameron.jack@anu.edu.au'
__status__ = 'Pre-release'


pos_args = ['db_path']
opt_args = ['sample']

script_info = {}
script_info['title'] = 'Summarises all DB info'
script_info['script_description'] = 'Creates a report of DB contents'
script_info['brief_info'] = 'Reports DB contents'
script_info['version'] = __version__
script_info['authors'] = __author__
script_info['args'] = Args(positional_args=pos_args, optional_args=opt_args)
script_info['required_options'] = script_info['args'].getReqCogentOpts()
script_info['optional_options'] = script_info['args'].getOptCogentOpts()

def main():
    rr = RunRecord('db_summary')
    rr.addCommands(sys.argv)
    args = script_info['args'].parse(window_title='DB Summary')
    session = make_session(args.db_path)
    sample_name = args.sample if args.sample else None

    chroms = get_chroms(session)
    species = get_species(session)

    if sample_name is None:
        total_samples_count = get_sample_counts(session)
        sample_names = get_all_sample_names(session)
        total_genes_count = get_gene_counts(session)
        total_exon_count = get_exon_counts(session)
        total_expr_count = get_expression_counts(session)
        total_diff_genes_count= get_diff_counts(session)
        total_target_genes_count = get_targetgene_counts(session)
        total_reffiles_count = get_reffile_counts(session)
    else:
        total_expr_count = get_expression_counts(session, sample_name)
        total_diff_genes_count= get_diff_counts(session, sample_name)
        total_target_genes_count = get_targetgene_counts(session, sample_name)
        reffiles_entries = get_reffile_entries(session, sample_name=sample_name)

    rr.addInfo('ChipPy DB name', args.db_path)
    rr.addInfo('Species name', species)
    rr.addInfo('Chroms list', chroms)
    if sample_name is None:
        rr.addInfo('Total # of sample entries', total_samples_count)
        rr.addInfo('Sample names', sample_names)
        rr.addInfo('Total # of gene entries', total_genes_count)
        rr.addInfo('Total # of exon entries', total_exon_count)
    rr.addInfo('Total # of absolute-scored gene entries', total_expr_count)
    rr.addInfo('Total # of differential gene entries', total_diff_genes_count)
    rr.addInfo('Total # of target gene entries', total_target_genes_count)
    if sample_name is None:
        rr.addInfo('Total # of reference files', total_reffiles_count)
    else:
        if len(reffiles_entries) > 0:
            rr.addInfo('Reference file name', reffiles_entries)
        else:
            rr.addError('Reference file name', 'Not Available')

    rr.display()

if __name__ == '__main__':
    main()

