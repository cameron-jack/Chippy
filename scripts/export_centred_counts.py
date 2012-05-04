"""dumps a RegionCollection compressed file"""
from __future__ import division

import os, sys, glob, warnings
warnings.filterwarnings('ignore', 'Not using MPI as mpi4py not found')
sys.path.extend(['..', '../src'])

import numpy

from optparse import make_option
from cogent.util.misc import parse_command_line_parameters

from chippy.core.count_tags import centred_counts_for_genes,\
            centred_diff_counts_for_genes,\
            centred_counts_external_genes
from chippy.core.collection import RegionCollection
from chippy.express import db_query
from chippy.express.db_schema import make_session
from chippy.draw.plottable import PlottableGroups
from chippy.ref.util import chroms
from chippy.util.run_record import RunRecord

__author__ = "Gavin Huttley, Cameron Jack"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley, Cameron Jack"
__credits__ = ["Gavin Huttley, Cameron Jack"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"
__version__ = '0.1'

# TODO: fix hard-wiring to Mouse
def get_collection(session, sample_name, expr_area, counts_dir, max_read_length,
        count_max_length, window_size, multitest_signif_val, filename,
        overwrite, sample_type, tab_delimited, run_record=None, test_run=False):

    if run_record is None:
        run_record = RunRecord()

    if not os.path.exists(filename) or overwrite:
        data_collection = None
        if sample_type == 'Expression data: absolute ranked':
            print "Collecting data for absolute expression experiment"
            data_collection, run_record = centred_counts_for_genes(session,
                    sample_name, expr_area, 'mouse', None, counts_dir, max_read_length,
                    count_max_length, window_size, run_record, test_run)
        
        elif sample_type == 'Expression data: difference in expression between samples':
            print "Collecting data for difference expression experiment"
            data_collection, run_record = centred_diff_counts_for_genes(session,
                    sample_name, expr_area, 'mouse', None, counts_dir, max_read_length,
                    count_max_length, window_size, multitest_signif_val,
                    run_record, test_run)
            
        else:
            print "Experiment type %s not supported" % sample_type
            
        if data_collection is not None:
            data_collection.writeToFile(filename, as_table=tab_delimited, compress_file=True)
        else:
            sys.stderr.write('No data_collection was returned!\n')
    else:
        print 'Existing output at %s' % filename
        data_collection = RegionCollection(filename=filename)
    
    return data_collection, run_record

if 'CHIPPY_DB' in os.environ:
    db_path = os.environ['CHIPPY_DB']
else:
    raise RuntimeError('You need to set an environment variable CHIPPY_DB '\
                       'that indicates where to find the database')

session = make_session( "sqlite:///%s" % db_path)

samples = db_query.get_samples(session)
if not samples:
    samples = [None]

script_info = {}
script_info['title'] = 'Saves feature centred counts'
script_info['script_description'] = 'Saves centred counts for TSS and '\
        'Exon-3prime, Intron-3prime or Exon 3&5-prime boundaries for a '\
        'given window size'
script_info['version'] = __version__
script_info['authors'] = __author__
script_info['output_description']= 'Generates a Pickle file or a gzipped '\
        'tab-delimited file that can be used for plotting of subsets of genes'

# options organisation

# essential sample specification
opt_sample = make_option('-c', '--sample', type='choice',
           help='Choose the expression study [default: %default]',
           choices=[str(s) for s in samples])

exp_absolute = 'Expression data: absolute ranked'
exp_diff = 'Expression data: difference in expression between samples'
external_genes ='External gene list'

opt_sample_type = make_option('-y', '--sample_type', type='choice',
        choices=[exp_absolute, exp_diff, external_genes],
        help='Select the type of data you want entered from %s' % \
                str([exp_absolute, exp_diff, external_genes]))

opt_expression_area = make_option('-e', '--expression_area', type='choice',
        choices=['TSS', 'Exon_3p', 'Intron_3p', 'Both_3p'], help='Expression '\
                'area options: TSS, Exon_3p, Intron-3p, Both-3p')

# essential source files
opt_counts_dir = make_option('-r', '--counts_dir',
    help='directory containing read counts. Can be a glob pattern for '\
    'multiple directories (e.g. for Lap1, Lap2 use Lap*)')
opt_save = make_option('-s', '--collection',
  help='path to save the plottable collection data '\
       +'(e.g. samplename-readsname-windowsize.gz)')

opt_overwrite = make_option('-f', '--overwrite',
             action='store_true', help='Ignore any saved files',
             default=False)

# optional counts generation
opt_read_length = make_option('-x', '--max_read_length', type='int',
        default=75, help='Maximum sequence read length [default: %default]')

opt_count_max_length = make_option('-k', '--count_max_length',
        action='store_false',
        help='Use maximum read length instead of mapped length',
        default=True)

opt_window = make_option('-w', '--window_size', type='int', default=1000,
        help='Region size around TSS [default: %default]')

opt_multitest_signif = make_option('-m', '--multitest_signif_val', type='int',
        help='Restrict plot to genes that pass multitest signficance,'\
        'valid values: 1, 0, -1', default=None)

#optional output
opt_tab_delimited = make_option('-d', '--tab_delimited', action='store_true',
        help='output to tab delimited format', default=False)

opt_test_run = make_option('-t', '--test_run', action='store_true',
        help="Test run, don't write output", default=False)

# adding into the main script_info dictionary required for correct processing
# via command-line or PyCogent.app
script_info['required_options'] = [opt_sample, opt_counts_dir, opt_save,
                                   opt_sample_type, opt_expression_area]

run_opts = [opt_overwrite, opt_tab_delimited, opt_test_run]
sampling_opts = [opt_read_length, opt_count_max_length, opt_window,
                 opt_multitest_signif]

script_info['optional_options'] = run_opts+sampling_opts

script_info['optional_options_groups'] = [('Run control', run_opts),
                                  ('Sampling', sampling_opts)
                                  ]

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    
    if opts.sample is None:
        raise RuntimeError('No samples available')
    
    sample_name = opts.sample.split(' : ')[0]
    print "Loading counts data for '%s'" % sample_name
    counts_dirs = opts.counts_dir
    dirname = os.path.dirname(counts_dirs)
    basename = os.path.basename(counts_dirs)
    counts_dirs = [os.path.join(dirname, p) for p in glob.glob1(dirname,
                                                    basename)]

    sample_type = opts.sample_type

    if (opts.multitest_signif_val is not None) and not (-1 <= opts.multitest_signif_val <= 1):
        raise RuntimeError('multitest_signif_val is not -1, 0, 1 or None. Halting execution.')

    if sample_type == exp_absolute or exp_diff:
        data_collection, rr = get_collection(session, sample_name,
                opts.expression_area, counts_dirs, opts.max_read_length,
                opts.count_max_length, opts.window_size,
                opts.multitest_signif_val, opts.collection, opts.overwrite,
                opts.sample_type, opts.tab_delimited, run_record=None,
                test_run=opts.test_run)

    else:
        print 'Other options not defined yet, choose from %s '\
              'or %s' % exp_absolute, exp_diff
    session.close()

    rr.display()

if __name__ == '__main__':
    main()

