#!/usr/bin/env python
import sys, os
from optparse import make_option

from cogent.util.misc import parse_command_line_parameters
from cogent import LoadTable

from chippy.express.db_query import get_samples
from chippy.express.db_populate import add_ensembl_gene_data, \
        add_expression_study, add_expression_diff_study, add_samples,\
        ReferenceFile, add_external_genes, add_samples
from chippy.express.db_schema import make_session
from chippy.parse.r_dump import SimpleRdumpToTable

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"
__version__ = '0.1'

if 'CHIPPY_DB' in os.environ:
    db_path = os.environ['CHIPPY_DB']
else:
    raise RuntimeError('You need to set an environment variable CHIPPY_DB '\
                       'that indicates where to find the database')

session = make_session( "sqlite:///%s" % db_path)

script_info = {}

script_info['title'] = 'Add expression study'
script_info['script_description'] = "Add an expression study from an R "\
                                    "export."
script_info['version'] = __version__

script_info['required_options'] = [
 make_option('-e',
         '--expression_data',
          help='Path to the expression data file. Must be tab delimited.'),
]


gene_id = make_option('-g','--gene_id_heading',
        default='ENSEMBL',
        help='Column containing the Ensembl gene stable ID')
probesets = make_option('-p','--probeset_heading',
        default='probeset',
        help='Column containing the probeset IDs')
exp_scores = make_option('-o','--expression_heading',
        default='exp',
        help='Column containing the expression scores')
identifiers = [gene_id, probesets, exp_scores]

_samples = get_samples(session)
if _samples:
    _samples = [str(s) for s in _samples]
else:
    _samples = [None]

existing_sample = make_option('-s','--sample', type='choice',
        choices=_samples,
        help='Select an existing or use field below to add new',
        )
default_new_sample = 'samplename : a description'
new_sample = make_option('-S','--new_sample',
        default=default_new_sample,
        help="Replace the text on the left and right of the ', "\
        "e.g. `S : S phase'")


exp_absolute = 'Expression data: absolute ranked'
exp_diff = 'Expression data: difference in expression between samples'
external_genes ='External gene list'

sample_type = make_option('-y', '--sample_type', type='choice',
        choices= [exp_absolute, exp_diff, external_genes],
            help='Select whether the type of data you want entered')


sample = [existing_sample, new_sample, sample_type]
related_file1 = make_option('--reffile1', default=None,
          help='Related file 1')
related_file2 = make_option('--reffile2', default=None,
          help='Related file 2')

related_files = [related_file1, related_file2]

script_info['optional_options'] = identifiers + sample + related_files

script_info['optional_options_groups'] = [
            ('Data file identifier columns (expression data)', identifiers),
            ('Sample', sample),
            ('Expression difference: related reference files', related_files)
            ]


def _get_name_description(value):
    """returns the name and description of a : separated sample"""
    new_sample = value.split(':')
    assert len(new_sample) == 2, "Only one ':' allowed!"
    name, description = new_sample
    name = name.strip()
    description = description.strip()
    return name, description

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    
    if opts.new_sample != default_new_sample:
        name, description = _get_name_description(opts.new_sample)
        rr = add_samples(session, [(name, description)])
    elif opts.sample is None:
        raise RuntimeError('No sample specified')
    else:
        name, description = _get_name_description(opts.sample)
        rr = None
    
    sample_type = opts.sample_type
    
    if sample_type in (exp_absolute, exp_diff):
        data = SimpleRdumpToTable(opts.expression_data,
                stable_id_label=opts.gene_id_heading,
                probeset_label=opts.probeset_heading,
                exp_label=opts.expression_heading, validate=True)
    else:
        raise NotImplementedError
    
    # TODO seems I should define constants for the headers and do a conversion
    # on reading in. This would eliminate all the optional arguments for
    # column labels
    if sample_type == exp_absolute:
        rr = add_expression_study(session, name, opts.expression_data, data,
                probeset_label=opts.probeset_heading, run_record=rr)
    elif sample_type == exp_diff:
        # diff between two files, check we got the related files
        assert opts.related_file1 is not None and \
            opts.related_file2 is not None, 'To enter differences in gene '\
            'expression you must specify the 2 files that contain the '\
            'absolute measures.'
        rr = add_expression_diff_study(session, opts.expression_data, data,
            opts.related_file1, opts.related_file2,
            probeset_label=opts.probeset_heading,
            ensembl_id_label=opts.gene_id_heading,
            expression_label=opts.expression_heading,
            prob_label='rawp', sig_label='sig', run_record=rr)
    elif sample_type == external_genes:
        raise NotImplementedError
    else:
        raise RuntimeError('Unknown sample type')
    
    rr.display()
    
    session.commit()

if __name__ == "__main__":
    main()

