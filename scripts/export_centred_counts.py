"""dumps a RegionCollection compressed file"""
from __future__ import division

import os, sys, glob, warnings
warnings.filterwarnings('ignore', 'Not using MPI as mpi4py not found')
sys.path.extend(['..', '../src'])

from chippy.core.count_tags import centred_counts_for_genes,\
            centred_diff_counts_for_genes
from chippy.core.collection import RegionCollection
from chippy.express import db_query
from chippy.express.util import sample_types
from chippy.util.run_record import RunRecord
from chippy.util import command_args

__author__ = "Gavin Huttley, Cameron Jack"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley, Cameron Jack"
__credits__ = ["Gavin Huttley, Cameron Jack"]
__license__ = "GPL"
__maintainer__ = "Cameron Jack"
__email__ = "cameron.jack@anu.edu.au"
__status__ = "Pre-release"
__version__ = '0.1'

def set_environment():
    """ create all command-line option groups and set script_info """
    script_info = {}
    script_info['title'] = 'Saves feature centred counts'
    script_info['script_description'] = 'Saves centred counts for TSS and '\
                                        'Exon-3prime, Intron-3prime or Exon 3&5-prime boundaries for a '\
                                        'given window size'
    script_info['version'] = __version__
    script_info['authors'] = __author__
    script_info['output_description']= 'Generates a Pickle file or a gzipped '\
                                       'tab-delimited file that can be used for plotting of subsets of genes'

    if 'CHIPPY_DB' in os.environ:
        db_path = os.environ['CHIPPY_DB']
    else:
        raise RuntimeError('You need to set an environment variable '
                'CHIPPY_DB that indicates where to find the database')

    req_args = ['sample', 'sample_type', 'expression_area',
            'counts_dir',  'collection']
    opt_args = ['overwrite', 'tab_delimited', 'max_read_length',
                'count_max_length', 'window_size', 'multitest_signif_val',
                'include_target', 'exclude_target', 'test_run']

    inputs = command_args.Args(required_args=req_args,
            optional_args=opt_args, db_path=db_path)

    return inputs.parsed_args, db_path, script_info

def get_collection(session, sample_name, expr_area, species, counts_dir, max_read_length,
        count_max_length, window_size, multitest_signif_val, filename, overwrite,
        sample_type, tab_delimited, include_target=None, exclude_target=None,
        run_record=None, test_run=False):

    if run_record is None:
        run_record = RunRecord()

    if not os.path.exists(filename) or overwrite:
        data_collection = None
        if sample_type == sample_types['exp_absolute']:
            print "Collecting data for absolute expression experiment"
            data_collection, run_record = centred_counts_for_genes(session,
                    sample_name, expr_area, species, None, counts_dir,
                    max_read_length, count_max_length, window_size,
                    include_target, exclude_target, run_record, test_run)
        
        elif sample_type == sample_types['exp_diff']:
            print "Collecting data for difference expression experiment"
            data_collection, run_record = centred_diff_counts_for_genes(
                    session, sample_name, expr_area, species, None,
                    counts_dir, max_read_length, count_max_length,
                    window_size, multitest_signif_val, include_target,
                    exclude_target, run_record, test_run)
            
        else:
            print "Experiment type %s not supported" % sample_type
            
        if data_collection is not None:
            data_collection.writeToFile(filename, as_table=tab_delimited,
                    compress_file=True)
        else:
            sys.stderr.write('No data_collection was returned!\n')
    else:
        print 'Existing output at %s' % filename
        data_collection = RegionCollection(filename=filename)
    
    return data_collection, run_record

def main():
    """ Returns a pickle of size window length containing chromatin mapping
     averages per base, one per gene, ranked by expression """
    args, db_path, script_info = set_environment()
        
    if args.sample is None:
        raise RuntimeError('No samples available')

    db_name = str(db_path).split('/')[-1]
    species = db_name.split('_')[-1].rstrip('.db')
    
    sample_name = args.sample.split(' : ')[0]
    print "Loading counts data for '%s'" % sample_name
    counts_dirs = args.counts_dir
    dirname = os.path.dirname(counts_dirs)
    basename = os.path.basename(counts_dirs)
    counts_dirs = [os.path.join(dirname, p) for p in glob.glob1(dirname,
                                                    basename)]
    sample_type = args.sample_type

    include_name = None
    exclude_name = None
    if args.include_target:
        include_name = args.include_target.split(' : ')[0]

    if args.exclude_target:
        exclude_name = args.exclude_target.split(' : ')[0]

    if (args.multitest_signif_val is not None) and not (-1 <= args.multitest_signif_val <= 1):
        raise RuntimeError('multitest_signif_val is not -1, 0, 1 or None. Halting execution.')

    session = db_query.make_session('sqlite:///' + str(db_path))
    if sample_type in [sample_types['exp_absolute'], sample_types['exp_diff']]:
        data_collection, rr = get_collection(session, sample_name,
                args.expression_area, species, counts_dirs, args.max_read_length,
                args.count_max_length, args.window_size,
                args.multitest_signif_val, args.collection, args.overwrite,
                args.sample_type, args.tab_delimited, include_name,
                exclude_name, run_record=None, test_run=args.test_run)

    else:
        print 'Other options not defined yet, choose from %s '\
                'or %s' % sample_types['exp_absolute'], \
                sample_types['exp_diff']

    session.close()
    rr.display()

if __name__ == '__main__':
    main()

