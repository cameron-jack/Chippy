#!/usr/bin/env python
"""Takes an R-dump file and populates DB - requires start_chippy_db has been run"""
import sys
sys.path.append('..')

import os
from cogent import LoadTable

from chippy.express import db_query
from chippy.express.db_populate import add_data
from chippy.express.util import sample_types
from chippy.parse.r_dump import SimpleRdumpToTable
from chippy.util import command_args
from chippy.util.run_record import RunRecord

__author__ = "Gavin Huttley, Cameron Jack"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley, Cameron Jack"
__credits__ = ["Gavin Huttley, Cameron Jack"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"
__version__ = '0.2'

def set_environment():
    """ create the DB session and run options """

    if 'CHIPPY_DB' in os.environ:
        db_path = os.environ['CHIPPY_DB']
    else:
        raise RuntimeError('You need to set an environment variable '
                           'CHIPPY_DB that indicates where to find the database')

    script_info = {}
    script_info['title'] = 'Add expression study'
    script_info['script_description'] = 'Add an expression study or a '\
            'differential expression study from a tab delimited file, or '\
            'load an ENSEMBL gene list.'

    script_info['version'] = __version__
    script_info['authors'] = __author__
    script_info['output_description']= 'None.'

    # Process command-line arguments
    req_args = ['expression_data', 'sample', 'new_sample',
                    'sample_type']
    opt_args = ['reffile1', 'reffile2', 'allow_probeset_many_gene',
                'gene_id_heading', 'probeset_heading', 'expression_heading',
                'test_run']
    inputs = command_args.Args(required_args=req_args,
            optional_args=opt_args, db_path=db_path)

    return inputs.parsed_args, db_path, script_info

def _get_name_description(value):
    """returns the name and description of a : separated sample"""
    new_sample = value.split(':')
    if new_sample[0].lower() == 'none':
        print 'No existing sample given'
        return 'None', 'None'
    assert len(new_sample) != 1, "Only one ':' allowed!"
    name, description = new_sample
    name = name.strip()
    description = description.strip()
    return name, description

def main():
    args, db_path, script_info = set_environment()
    session = db_query.make_session('sqlite:///' + db_path)
    rr = RunRecord()

    if args.new_sample.count(':') == 1:
        name, description = _get_name_description(args.new_sample)
    elif args.sample is None:
        raise RuntimeError('No sample specified')
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
        ensembl_id_label=args.gene_id_heading,
        expr_heading=args.expression_heading,
        sample_type=args.sample_type, reffile1=args.reffile1,
        reffile2=args.reffile2, rr=RunRecord())

    rr.addInfo('add_expression_db', name+' added to DB', success)
    rr.display()

if __name__ == "__main__":
    main()
