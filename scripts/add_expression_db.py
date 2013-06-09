#!/usr/bin/env python
"""Takes an R-dump file and populates DB - requires start_chippy_db has been run"""
import sys
sys.path.append('..')

from cogent import LoadTable

from chippy.express import db_query
from chippy.express.db_populate import add_data
from chippy.express.util import sample_types
from chippy.parse.r_dump import simpleRdumpToTable
from chippy.util.command_args import Args
from chippy.util.run_record import RunRecord

__author__ = 'Gavin Huttley, Cameron Jack'
__copyright__ = 'Copyright 2011-2013, Anuj Pahwa, Gavin Huttley, Cameron Jack'
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
req_args = ['expression_data', 'new_sample', 'sample_type']
opt_args = ['reffile1', 'reffile2', 'allow_probeset_many_gene',
        'gene_id_heading', 'probeset_heading', 'expression_heading',
        'test_run']
pos_args = ['db_path']

script_info['args'] = Args(required_args=req_args, optional_args=opt_args,
    positional_args=pos_args)
script_info['required_options'] = script_info['args'].req_cogent_opts
script_info['optional_options'] = script_info['args'].opt_cogent_opts

def main():
    rr = RunRecord('add_expression_db')
    rr.addCommands(sys.argv)

    args = script_info['args'].parse()
    session = db_query.make_session(args.db_path)

    if args.new_sample is None:
        raise RuntimeError('No sample specified')
    elif args.new_sample.count(':') == 1:
        name, description = args.new_sample.split(':')
        name = name.strip()
        description = description.strip()
    else:
        raise RuntimeError("Please provide 'Name : Description'")

    ref_file = args.expression_data

    # Check that ReferenceFile is unique
    if db_query.get_reffile_counts(session, reffile_name=ref_file) == 1:
        rr.addCritical('ReferenceFile already loaded', ref_file)
        rr.display()
        sys.exit(1)

    if sample_types[args.sample_type] in\
            (sample_types['exp_absolute'], sample_types['exp_diff']):
        expr_table = simpleRdumpToTable(args.expression_data,
                stable_id_label=args.gene_id_heading,
                probeset_label=args.probeset_heading,
                exp_label=args.expression_heading, validate=True)
    else:
        expr_table = LoadTable(args.expression_data, sep='\t')

    success = add_data(session, name, description,
            args.expression_data, expr_table,
            probeset_heading=args.probeset_heading,
            gene_id_heading=args.gene_id_heading,
            expr_heading=args.expression_heading,
            sample_type=args.sample_type, reffile1=args.reffile1,
            reffile2=args.reffile2)

    rr.addInfo(name+' added to DB', success)
    rr.display()

if __name__ == "__main__":
    main()
