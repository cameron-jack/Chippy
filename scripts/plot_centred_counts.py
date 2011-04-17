from __future__ import division

import os, sys
sys.path.extend(['../../src'])

import numpy

from optparse import make_option
from cogent.util.misc import parse_command_line_parameters

from chippy.core.count_tags import centred_counts_for_genes,\
            centred_counts_external_genes
from chippy.core.collection import RegionCollection
from chippy.express import db_query
from chippy.draw.plottable import PlottableGroups
from chippy.ref.util import chroms
from chippy.util.util import create_path, make_cl_command, just_filename, \
                    dirname_or_default

__author__ = 'Gavin Huttley'
__copyright__ = 'Copyright 2011, Gavin Huttley'
__credits__ = ['Gavin Huttley']
__license__ = 'GPL'
__maintainer__ = 'Gavin Huttley'
__email__ = 'Gavin.Huttley@anu.edu.au'
__status__ = 'alpha'
__version__ = '0.1'

def make_sample_choices(session):
    """returns the available choices for external gene samples"""
    samples = ['%s : %s' % (s.name, s.description)
        for s in db_query.get_external_sample(session)]
    samples.insert(0, None)
    return samples

def get_sample_name(sample):
    """returns sample name from a 'sample : description' string"""
    if str(sample) != 'None':
        sample = sample.split(':')[0].strip()
    else:
        sample = None
    return sample

if 'CHIPPY_DB' not in os.environ:
    # TODO remove this hardcoding when PyCogent.app allows specifying
    # environment variables. This should be an assert statement, and fail
    # if not present in the environment.
    db_path = "/Users/gavin/Desktop/Pending/lap1/chippy.db"
else:
    db_path = os.environ['CHIPPY_DB']

ensembl_release='58'
session = db_query.make_session('sqlite:///%s' % db_path)
samples = db_query.get_external_sample(session)

script_info = {}
script_info['title'] = 'Plot read counts heat-mapped by gene expression'
script_info['script_description'] = "Takes read counts that are centred on"\
    " on gene TSS, sorted from high to low gene expression and makes a"\
    " heat-mapped plot."
script_info['version'] = __version__
script_info['authors'] = __author__
script_info['output_description']= "Generates either a single pdf figure or"\
    " a series of pdfs that can be merged into a movie."

# alternate option organisation

# essential source files
opt_collection = make_option('-s', '--collection',
  help='path to the plottable data '\
       +'(e.g. samplename-readsname-windowsize.gz)')

# chrom choice
opt_chroms = make_option('-C', '--chrom', type='choice', default='All',
               help='Choose a chromosome [default: %default]',
               choices=('All',)+chroms['mouse'])

# or external sample (gene) choice
opt_extern = make_option('-E', '--external_sample', type='choice', 
            default=None,
            choices=make_sample_choices(session),
            help='External sample')

# essential plotting information
opt_grp_size = make_option('-g', '--group_size', type='choice', default='All',
   choices=['All', '50', '100'],
   help='Number of genes to group to estimate statistic [default: %default]')

# optional plotting information
opt_cutoff = make_option('-k', '--cutoff', type='float', default = 3,
             help='Stdev limit [default: %default]')
opt_fig_height = make_option('-H', '--fig_height', type='float', default=2.5*6,
   help='Figure height (cm) [default: %default]')

opt_fig_width = make_option('-W', '--fig_width', type='float', default=2.5*12,
   help='Figure width (cm) [default: %default]')
opt_bgcolor = make_option('-b', '--bgcolor', type='choice', default='black',
               help='Plot background color [default: %default]',
               choices=['black', 'white'])
opt_yrange = make_option('-y', '--ylim', default=None,
       help='minimum-maximum yaxis values (e.g. 0-3.5)')
opts_xgrid_locate = make_option('--xgrid_lines', type='float', default = 100,
                 help='major grid-line spacing on x-axis [default: %default]')
opts_ygrid_locate = make_option('--ygrid_lines', type='float', default = 0.5,
                 help='major grid-line spacing on y-axis [default: %default]')
opts_xlabel_interval = make_option('--xlabel_interval', type='int',
        default = 2,
        help='number of blank ticks between labels [default: %default]')
opts_ylabel_interval = make_option('--ylabel_interval', type='int',
    default = 2,
    help='number of blank ticks between labels [default: %default]')
opts_xlabel_font = make_option('--xfontsize', type='int',
    default = 12,
    help='font size for x label [default: %default]')
opts_ylabel_font = make_option('--yfontsize', type='int',
    default = 12,
    help='font size for y label [default: %default]')
opts_vline_style = make_option('--vline_style', type='choice',
    default = '-.', choices=['-.', '-', '.'],
    help='line style for centred vertical line [default: %default]')
opts_vline_width = make_option('--vline_width', type='int',
    default = 2, 
    help='line width for centred vertical line [default: %default]')
opts_ylabel = make_option('--ylabel',
    default = 'Mean counts', help='Label for the y-axis [default: %default]')
opts_xlabel = make_option('--xlabel',
    default = 'Position relative to TSS',
    help='Label for the x-axis [default: %default]')
opts_title = make_option('--title', help='Plot title [default: %default]')
opts_alpha = make_option('--line_alpha', type='float', default=1.0,
                 help='Opacity of lines [default: %default]')
opts_plot_filename = make_option('--plot_filename',
    default = None,
    help='Name of final plot file (must end with .pdf) [default: %default]')

# 
opts_plotseries = make_option('-p', '--plot_series',
                 action='store_true', default=False,
     help='Plot series of figures. A directory called plot_filename-series'\
          +' will be created. Requires plot_filename be defined.')

opt_txt_coords = make_option('--text_coords', default=None,
       help='x, y coordinates of series text (e.g. 600,3.0)')

# 
opt_test_run = make_option('-t', '--test_run',
             action='store_true', help="Test run, don't write output",
             default=False)

script_info['required_options'] = [opt_collection, opt_yrange]

run_opts = [opt_test_run]
sampling_opts = [opt_grp_size, opt_extern, opt_chroms, opt_cutoff]
save_opts = [opts_plot_filename]
series_opts = [opts_plotseries, opt_txt_coords]
plot_labels = [opts_title, opts_ylabel, opts_xlabel]
plot_dims = [opt_fig_height, opt_fig_width, opts_xgrid_locate,
            opts_ygrid_locate, opts_xlabel_interval, opts_ylabel_interval]
plot_colors = [opt_bgcolor, opts_alpha, opts_vline_style, opts_vline_width,
        opts_xlabel_font, opts_ylabel_font]

script_info['optional_options'] = run_opts+sampling_opts+save_opts+\
        series_opts+plot_labels+plot_dims+plot_colors

script_info['optional_options_groups'] = [('Run control', run_opts),
                                  ('Sampling', sampling_opts),
                                  ('Saving', save_opts),
                                  ('Plot series', series_opts),
                                  ('Plot text', plot_labels),
                                  ('Plot colours', plot_colors),
                                  ('Plot dimensions', plot_dims)
                                  ]


def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    
    ylim = map(float, opts.ylim.strip().split('-'))
    print 'Loading counts data'
    data_collection = RegionCollection(filename=opts.collection)
    total_gene = data_collection.ranks.max()
    
    external_sample = get_sample_name(opts.external_sample)
    ensembl_release = data_collection.info['args']['ensembl_release']
    stable_ids = None
    if external_sample is not None:
        genes = db_query.get_external_genes(session, ensembl_release,
                                            external_sample)
        stable_ids = [g.ensembl_id for g in genes]
    elif opts.chrom != 'All':
        genes = db_query.get_genes(session, ensembl_release, opts.chrom)
        stable_ids = [g.ensembl_id for g in genes]
    
    if stable_ids is not None:
        data_collection = data_collection.filteredByLabel(stable_ids)
        session.close()
    
    data_collection.ranks /= total_gene
    print '\tNumber of sampled genes: %d' % data_collection.N
    
    window_size = data_collection.info['args']['window_size']
    # exclude genes with a count > opts.cutoff std above mean
    # x here will be a normalised statistic -- but only if the genes are not
    # from an external_sample
    if external_sample is None:
        data_collection = data_collection.\
                                    filteredNormalised(cutoff=opts.cutoff)
    
    print 'Getting grouped data'
    # if we have a plot series, we need to create a directory to dump the
    # files into
    if opts.plot_series:
        save_dir = dirname_or_default(opts.plot_filename)
        basename = os.path.basename(opts.plot_filename)
        plot_filename = os.path.join(save_dir, basename)
        
        plot_series_dir = os.path.join(save_dir,
                        '%s-series' % basename[:basename.rfind('.')])
        create_path(plot_series_dir)
        labels = []
        filename_series = []
    else:
        labels = None
        filename_series = None
        series_labels = None
        label_coords = None
    
    counts = []
    ranks = []
    group_size = [opts.group_size, data_collection.N][opts.group_size=='All']
    group_size = int(group_size)
    for index, (c,r,l) in enumerate(data_collection.iterTransformedGroups(
                        group_size=group_size)):
        counts.append(c)
        ranks.append(r)
        if opts.plot_series:
            labels.append('Group %d' % index)
    
    num_groups = len(counts)
    print '\tNumber of groups: %d' % num_groups
    # reverse the counts and colour series so low color goes first
    counts = list(reversed(counts))
    ranks = list(reversed(ranks))
    if opts.plot_series:
        label_coords = map(float, opts.text_coords.split(','))
        series_labels = list(reversed(labels))
        series_template = 'plot-%%.%sd.pdf' % len(str(len(counts)))
        filename_series = [os.path.join(plot_series_dir, series_template % i)
                            for i in range(len(series_labels))]
    
    print 'Starting to plot'
    if opts.bgcolor == 'black':
        grid={'color': 'w'}
        bgcolor='0.1'
        vline_color='w'
    else:
        grid={'color': 'k'}
        vline_color='k'
        bgcolor='0.9'
    
    vline = dict(x=0, linewidth=opts.vline_width,
                   linestyle=opts.vline_style, color=vline_color)
    
    plot = PlottableGroups(height=opts.fig_height/2.5,
        width=opts.fig_width/2.5,
        bgcolor=bgcolor, grid=grid,
        ylim=ylim, xlim=(-window_size, window_size),
        xtick_space=opts.xgrid_lines, ytick_space=opts.ygrid_lines,
        xtick_interval=opts.xlabel_interval,
        ytick_interval=opts.ylabel_interval,
        xlabel_fontsize=opts.xfontsize, ylabel_fontsize=opts.yfontsize,
        vline=vline,
        ioff=True)
    
    x = numpy.arange(-window_size, window_size)
    plot(x, y_series=counts, color_series=ranks, series_labels=series_labels,
        filename_series=filename_series, label_coords=label_coords,
        alpha=opts.line_alpha, xlabel='Position relative to TSS',
        ylabel='Mean counts', title = opts.title)
    plot.savefig(opts.plot_filename)
    plot.show()
    

if __name__ == '__main__':
    main()

