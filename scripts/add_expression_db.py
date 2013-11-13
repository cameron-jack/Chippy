#!/usr/bin/env python
"""Takes an R-dump file and populates DB - requires start_chippy_db has been run"""
import sys
sys.path.append('..')

from cogent import LoadTable

from chippy.express import db_query
from chippy.express.db_populate import add_data
from chippy.express.util import sample_types
from chippy.parse.expr_data import gene_expr_to_table, gene_expr_diff_to_table
from chippy.util.command_args import Args
from chippy.util.run_record import RunRecord

__author__ = 'Gavin Huttley, Cameron Jack'
__copyright__ = 'Copyright 2011-2013, Gavin Huttley, Cameron Jack, Anuj Pahwa'
__credits__ = ['Gavin Huttley', 'Cameron Jack']
__license__ = 'GPL'
__maintainer__ = 'Cameron Jack'
__email__ = 'cameron.jack@anu.edu.au'
__status__ = 'Pre-release'
__version__ = '0.2'

script_info = {}
script_info['title'] = 'Add expression study'
script_info['script_description'] = 'Add an expression study or a '+\
        'differential expression study from a tab delimited file, or '+\
        'load an ENSEMBL gene list as target genes.'
script_info['brief_description'] = 'Adds an expression study to a '+\
        'ChippyDB'
script_info['version'] = __version__
script_info['authors'] = __author__
script_info['output_description']= 'None.'

# Process command-line arguments
req_args = ['expression_data', 'name', 'description', 'sample_type']
opt_args = ['reffile1', 'reffile2', 'allow_probeset_many_gene',
        'gene_id_heading', 'probeset_heading', 'expression_heading',
        'significance_heading', 'p_value_heading', 'sep']
pos_args = ['db_path']

script_info['args'] = Args(required_args=req_args, optional_args=opt_args,
    positional_args=pos_args)
script_info['required_options'] = script_info['args'].req_cogent_opts
script_info['optional_options'] = script_info['args'].opt_cogent_opts

def main():
    rr = RunRecord('add_expression_db')
    rr.addCommands(sys.argv)

    args = script_info['args'].parse(use_scrollbars=True,
            use_save_load_button=True, window_title='Add Expression to DB')
    session = db_query.make_session(args.db_path)

    name = args.name
    description = args.description
    ref_file = args.expression_data
    sample_type=args.sample_type

    # Check that Sample and Reference File are both unique
    if name in db_query.get_sample_entries(session):
        rr.dieOnCritical('Sample name already exists', name)
    if ref_file in db_query.get_reffile_entries(session, reffile_name=ref_file):
        rr.dieOnCritical('ReferenceFile already loaded', ref_file)

    if sample_types[sample_type] == sample_types['abs_expr']:
        expr_table = gene_expr_to_table(args.expression_data,
                stable_id_label=args.gene_id_heading,
                probeset_label=args.probeset_heading,
                exp_label=args.expression_heading,
                allow_probeset_many_gene=args.allow_probeset_many_gene,
                validate=True, sep=args.sep)

    elif sample_types[sample_type] == sample_types['diff_expr']:
        # validation breaks with some of Rohan's diff files
        # he's included all probesets but only the mean score, once.
        expr_table = gene_expr_diff_to_table(args.expression_data,
                stable_id_label=args.gene_id_heading,
                probeset_label=args.probeset_heading,
                exp_label=args.expression_heading,
                sig_label=args.significance_heading,
                pval_label=args.p_value_heading,
                allow_probeset_many_gene=args.allow_probeset_many_gene,
                validate=False, sep=args.sep)
    elif sample_types[sample_type] == sample_types['target_genes']:
        expr_table = LoadTable(args.expression_data, sep=args.sep)

    else:
        rr.dieOnCritical('Unknown sample type', args.sample_type)

    success = add_data(session, name, description,
            args.expression_data, expr_table,
            sample_type=args.sample_type, reffile1=args.reffile1,
            reffile2=args.reffile2)

    rr.addInfo(name+' added to DB', success)
    rr.display()

if __name__ == "__main__":
    main()
