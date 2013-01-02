""" Converts legacy pickle ChIP-Seq data to BED format"""
from __future__ import division

import os, sys, glob, warnings
warnings.filterwarnings('ignore', 'Not using MPI as mpi4py not found')
sys.path.extend(['..', '../src'])

from optparse import make_option
from cogent.util.misc import parse_command_line_parameters
from cogent.util.progress_display import display_wrap
from cogent import LoadTable

from chippy.util.definition import NULL_STRAND, PLUS_STRAND, MINUS_STRAND

__author__ = "Cameron Jack, Gavin Huttley"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley, Cameron Jack"
__credits__ = ["Cameron Jack, Gavin Huttley"]
__license__ = "GPL"
__maintainer__ = "Cameron Jack"
__email__ = "cameron.jack@anu.edu.au"
__status__ = "Pre-release"
__version__ = '0.1'

script_info = {}
script_info['title'] = 'Converts all counts to BED format'
script_info['script_description'] = 'Converts all counts created by older pipeline to BED format'
script_info['version'] = __version__
script_info['authors'] = __author__
script_info['output_description']= 'BED6 file format. See http://genome.ucsc.edu/FAQ/FAQformat.html#format1'

# options organisation

# essential sample specification

# essential source files
opt_counts_dir = make_option('-r', '--counts_dir',
        help='directory containing read counts. Can be a glob pattern for '\
        'multiple directories (e.g. for Lap1, Lap2 use Lap*)')

opt_save = make_option('-s', '--save_path',
        help="path to save the output BED file "\
        +"(e.g. blah//samplename.bed)")

opt_overwrite = make_option('-f', '--force_overwrite',
        action='store_true', help='Ignore any saved files',
        default=False)

# optional counts generation
opt_read_length = make_option('-x', '--max_read_length', type='int',
        default=100, help='Maximum sequence read length [default: %default]')

opt_count_max_length = make_option('-k', '--count_max_length',
        action='store_false',
        help="Use maximum read length instead of mapped length",
        default=True)

opt_prefix = make_option('-p', '--prefix_chrom_number', type='string',
        default='chr', help='string to prefix chromosome number in BED '\
        +'[default: %default]')

opt_feature = make_option('--feature_name', type='string',
        help='string describing the mapped feature e.g. H2A.Z ')

opt_test_run = make_option('-t', '--test_run',
        action='store_true', help="Test run, don't write output",
        default=False)


# adding into the main script_info dictionary required for correct processing
# via command-line or PyCogent.app
script_info['required_options'] = [opt_counts_dir, opt_save, opt_feature]

run_opts = [opt_overwrite, opt_test_run]
sampling_opts = [opt_read_length, opt_count_max_length]

script_info['optional_options'] = run_opts+sampling_opts

script_info['optional_options_groups'] = [('Run control', run_opts),
        ('Sampling', sampling_opts)
]

@display_wrap
def make_bed_entries(mapped_read_path, chrom_number, feature_name, output_bed_file, max_read_length=None,
                       count_max_length=False, strand=NULL_STRAND, sep='\t',
                       is_sorted=True, ui=None):
    """ translates a numpy array of mapped chromosome positions into a BED file.

    Arguments:
        - mapped_read_path: path to table containing read coordinates, frequency data
        - output_write_path: path to BED-6 format output file
        - max_read_length: maximum length of a read length
        - count_max_length: if max_read_length provided, all mapped seqs set
          to this length
        - strand: only reads from specified strand are added. Default is both.
        - sep: the delimiter in the read coordinates file
        - is_sorted: whether the read file is already sorted
    """

    data = LoadTable(mapped_read_path, sep=sep)
    assert list(data.Header) == ['start', 'length', 'strand', 'freq'],\
    "mapped read Table header doesn't match expected"

    if not is_sorted:
        data = data.sorted(columns='start')

    if count_max_length:
        assert max_read_length, 'must specify max_read_length to use'\
                                ' count_max_length'
    data = data.array.astype(int)
    total_data = data.shape[0]

    chrom = "chr%s" % (chrom_number)
    score = 30
    name = feature_name
    for i, row in enumerate(data):
        if i % 10 == 0:
            ui.display('Converting mapped locations [%d / %d]' % (i, total_data), i / total_data)
        start = row[0] # move from 1-based BWA calls to 0-based BED
        end = start + row[1] # end is 1 beyond actual mapped end, so length = end - start
        strand = ['+','-'][(row[2] > 0) == 0] # 1,-1 converted to '+','-'
        #print "\n Row values: %d\t%d\t%d\t%d\n" % (row[0], row[1], row[2], row[3])
        for counts in range(row[3]):
            bed_string = '%s\t%s\t%s\t%s\t%s\t%s\n' % \
                         (chrom, start, end, name, score, strand)
            #print bed_string
            output_bed_file.write(bed_string)

    print "Converted %d mapped locations" % (total_data)

def main():
    option_parser, opts, args =\
    parse_command_line_parameters(**script_info)

    if opts.counts_dir is None:
        raise RuntimeError('No data available')

    if os.path.exists(opts.save_path) and opts.force_overwrite is False:
        raise RuntimeError ('Output BED file exists and overwriting not enabled. Exiting.\n')

    output_bed_file = open(opts.save_path, 'w')

    dir_list = os.listdir(opts.counts_dir)

    for chrom_file in dir_list:
        if chrom_file[-7:] == '.txt.gz':
            chrom_number = (chrom_file.lstrip('chr').rstrip('.txt.gz'))
            print 'Chromosome file name: %s, chromosome number: %s\n' % (chrom_file, chrom_number)

            read_file = os.path.join(opts.counts_dir, chrom_file)

            make_bed_entries(read_file, chrom_number, opts.feature_name, output_bed_file,
                    opts.max_read_length, opts.count_max_length, chrom_number)
        else:
            print "\nSkipping file: %s as not data file.\n" % (chrom_file)

    output_bed_file.close()


if __name__ == '__main__':
    main()