#!/usr/bin/env python
"""Takes an R-dump file and populates DB - requires start_chippy_db has been run"""
import sys
sys.path.append('..')

from cogent import LoadTable

from chippy.express import db_query
from chippy.express.db_populate import add_data
from chippy.express.util import sample_types
from chippy.parse.r_dump import SimpleRdumpToTable
from chippy.util.command_args import Args
from chippy.util.run_record import RunRecord

__author__ = "Gavin Huttley, Cameron Jack"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley, Cameron Jack"
__credits__ = ["Gavin Huttley, Cameron Jack"]
__license__ = "GPL"
__maintainer__ = "Cameron Jack"
__email__ = "cameron.jack@anu.edu.au"
__status__ = "Pre-release"
__version__ = '0.2'

def set_environment():
    """ create the DB session and run options """
    script_info = {}
    script_info['title'] = 'Add expression study'
    script_info['script_description'] = 'Add an expression study or a '\
            'differential expression study from a tab delimited file, or '\
            'load an ENSEMBL gene list as target genes.'

    script_info['version'] = __version__
    script_info['authors'] = __author__
    script_info['output_description']= 'None.'

    # Process command-line arguments
    req_args = ['expression_data', 'new_sample',
            'sample_type']
    opt_args = ['reffile1', 'reffile2', 'allow_probeset_many_gene',
            'gene_id_heading', 'probeset_heading', 'expression_heading',
            'test_run']
    pos_args = ['db_path']
    inputs = Args(required_args=req_args, optional_args=opt_args,
            positional_args=pos_args)

    return inputs.parsed_args, script_info

def main():
    args, script_info = set_environment()
    session = db_query.make_session('sqlite:///' + args.db_path)
    rr = RunRecord()

    if args.sample is None:
        raise RuntimeError('No sample specified')
    elif args.new_sample.count(':') == 1:
        name, description = args.new_sample.split(':')
        name = name.strip()
        description = description.strip()
    else:
        raise RuntimeError("Please provide 'Name : Description'")

    if sample_types[args.sample_type] in\
            (sample_types['exp_absolute'], sample_types['exp_diff']):
        expr_table, rr = SimpleRdumpToTable(args.expression_data,
                stable_id_label=args.gene_id_heading,
                probeset_label=args.probeset_heading,
                exp_label=args.expression_heading, validate=True,
                run_record=rr)
    else:
        expr_table = LoadTable(args.expression_data, sep='\t')

    success, rr = add_data(session, name, description,
        args.expression_data, expr_table,
        probeset_heading=args.probeset_heading,
        gene_id_heading=args.gene_id_heading,
        expr_heading=args.expression_heading,
        sample_type=args.sample_type, reffile1=args.reffile1,
        reffile2=args.reffile2, rr=RunRecord())

    rr.addInfo('add_expression_db', name+' added to DB', success)
    rr.display()

if __name__ == "__main__":
    main()
