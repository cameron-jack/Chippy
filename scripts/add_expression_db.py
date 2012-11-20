#!/usr/bin/env python
"""Takes an R-dump file and populates DB - requires start_chippy_db has been run"""
import sys
sys.path.append('..')

import os
from cogent import LoadTable

from chippy.express import db_query
from chippy.express.db_populate import add_ensembl_gene_data, \
        add_expression_study, add_expression_diff_study, add_samples,\
        ReferenceFile, add_target_genes, add_samples
from chippy.express.db_schema import make_session
from chippy.parse.r_dump import SimpleRdumpToTable
from chippy.util import command_args

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

    ### All inputs are divided into logical groupings
    inputs = command_args.Args()
    # sample and new_sample are mutually exclusive
    inputs.add_args(['expression_data', 'sample', 'new_sample',
            'sample_type'], db_path=db_path, required=True)
    inputs.add_args(['reffile1', 'reffile2', 'allow_probeset_many_gene',
            'gene_id_heading', 'probeset_heading', 'expression_heading',
            'test_run'])
    args = inputs.parser.parse_args()

    return args, db_path, script_info

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

    if args.new_sample.count(':') == 1:
        name, description = _get_name_description(args.new_sample)
        successes, rr = add_samples(session, [(name, description)])
        # We try to added a single sample and fail
        if len(successes) == 1 and successes[0] == False:
            rr.display()
            sys.exit(-1)

    elif args.sample is None:
        raise RuntimeError('No sample specified')
    else:
        name, description = _get_name_description(args.sample)
        rr = None

    sample_type = args.sample_type

    exp_absolute = 'Expression data: absolute ranked'
    exp_diff = 'Expression data: difference in expression between samples'
    target_genes = 'Target gene list'

    if sample_type in (exp_absolute, exp_diff):
        data, rr = SimpleRdumpToTable(args.expression_data,
                stable_id_label=args.gene_id_heading,
                probeset_label=args.probeset_heading,
                exp_label=args.expression_heading, validate=True,
                run_record=rr)
    else:
        data = LoadTable(args.expression_data, sep='\t')

    if sample_type == exp_absolute:
        rr = add_expression_study(session, name, args.expression_data, data,
            probeset_label=args.probeset_heading,
            ensembl_id_label=args.gene_id_heading, run_record=rr)
    elif sample_type == exp_diff:
        # diff between two files, check we got the related files
        assert args.reffile1 is not None and args.reffile2 is not None,\
            'To enter differences in gene expression you must specify the 2'\
            'files that contain the absolute measures.'
        rr = add_expression_diff_study(session, name,
            args.expression_data, data, args.reffile1, args.reffile2,
            probeset_label=args.probeset_heading,
            ensembl_id_label=args.gene_id_heading,
            expression_label=args.expression_heading,
            prob_label='rawp', sig_label='sig', run_record=rr)
    elif sample_type == target_genes:
        rr = add_target_genes(session, name, args.expression_data, data,
                    ensembl_id_label=args.gene_id_heading, run_record=rr)
    else:
        raise RuntimeError('Unknown sample type')
    
    rr.display()
    
    session.commit()

if __name__ == "__main__":
    main()

