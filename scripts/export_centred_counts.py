"""dumps a RegionCollection compressed file"""
from __future__ import division

import os, sys, glob
sys.path.extend(['..', '../src'])

import numpy

from optparse import make_option
from cogent.util.misc import parse_command_line_parameters

from chippy.core.count_tags import centred_counts_for_genes,\
            centred_counts_external_genes
from chippy.core.collection import RegionCollection
from chippy.express import db_query
from chippy.draw.plottable import PlottableGroups
from chippy.ref.util import chroms

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"
__version__ = '0.1'

def get_collection(session, sample_name, ensembl_release, counts_dir,
    max_read_length, count_max_length, window_size, filename, overwrite, test_run):
    if not os.path.exists(filename) or overwrite:
        data_collection = centred_counts_for_genes(session, sample_name,
                'mouse', None, counts_dir, ensembl_release, max_read_length,
                count_max_length, window_size, test_run)
        if data_collection is not None:
            data_collection.writeToFile(filename)
        else:
            sys.stderr.write('No data_collection was returned!\n')
    else:
        data_collection = RegionCollection(filename=filename)
    
    return data_collection

if 'CHIPPY_DB' in os.environ:
    db_path = os.environ['CHIPPY_DB']
else:
    raise RuntimeError('You need to set an environment variable CHIPPY_DB '\
                       'that indicates where to find the database')

# TODO remove this hardcoding to an ENSEMBL version
ensembl_release='58'
session = db_query.make_session('sqlite:///%s' % db_path)
samples = db_query.get_samples(session)
if not samples:
    samples = [None]

script_info = {}
script_info['title'] = 'Saves TSS centred counts'
script_info['script_description'] = 'Saves centred counts from a sample for a given window size'
script_info['version'] = __version__
script_info['authors'] = __author__
script_info['output_description']= 'Generates either a compressed file that can be used for plotting of subsets of genes'

# options organisation

# essential sample specification
opt_sample = make_option('-c', '--sample', type='choice',
           help='Choose the expression study [default: %default]',
           choices=[str(s) for s in samples])

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
             help="Use maximum read length instead of mapped length",
             default=True)

opt_window = make_option('-w', '--window_size', type='int', default=1000,
                help='Region size around TSS [default: %default]')

opt_test_run = make_option('-t', '--test_run',
             action='store_true', help="Test run, don't write output",
             default=False)

# adding into the main script_info dictionary required for correct processing
# via command-line or PyCogent.app
script_info['required_options'] = [opt_sample, opt_counts_dir, opt_save]

run_opts = [opt_overwrite, opt_test_run]
sampling_opts = [opt_read_length, opt_count_max_length, opt_window]

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
    print 'Loading counts data'
    counts_dirs = opts.counts_dir
    dirname = os.path.dirname(counts_dirs)
    basename = os.path.basename(counts_dirs)
    counts_dirs = [os.path.join(dirname, p) for p in glob.glob1(dirname,
                                                    basename)]
    data_collection = get_collection(session, sample_name, ensembl_release,
        counts_dirs , opts.max_read_length, opts.count_max_length, 
        opts.window_size, opts.collection, opts.overwrite, opts.test_run)
    session.close()

if __name__ == '__main__':
    main()

