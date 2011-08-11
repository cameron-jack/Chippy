from __future__ import division

import os, sys
sys.path.extend(['../../src'])

import numpy

from optparse import make_option
from cogent.util.misc import parse_command_line_parameters
from cogent.maths.stats.jackknife import JackknifeStats

from chippy.core.count_tags import centred_counts_for_genes,\
            centred_counts_external_genes
from chippy.core.collection import RegionCollection, column_sum, column_mean
from chippy.express import db_query
from chippy.draw.plottable import PlottableSingle
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

if 'CHIPPY_DB' in os.environ:
    db_path = os.environ['CHIPPY_DB']
else:
    raise RuntimeError('You need to set an environment variable CHIPPY_DB '\
                       'that indicates where to find the database')

def stat_maker(func, coll):
    def calc_stat(coords):
        subset_data = coll.take(coords)
        return func(subset_data)
    return calc_stat

def summed(data):
    freqs = data.asfreqs()
    c, r = freqs.transformed(counts_func=column_sum)
    return c

def averaged(data):
    c, r = data.transformed(counts_func=column_mean)
    return c

session = db_query.make_session('sqlite:///%s' % db_path)
samples = db_query.get_external_sample(session)

script_info = {}
script_info['title'] = 'Compare read counts between histone variants'
script_info['script_description'] = "Takes read counts that are centred on"\
    " on gene TSS, that exist in two separate mapped read samples."
script_info['version'] = __version__
script_info['authors'] = __author__
script_info['output_description']= "Generates a single pdf figure."

# alternate option organisation

# essential source files
opt_collection1 = make_option('-1', '--collection1',
  help='path to the plottable data from sample 1'\
       +'(e.g. samplename-readsname-windowsize.gz)')
opt_collection2 = make_option('-2', '--collection2',
  help='path to the plottable data from sample 2'\
       +'(e.g. samplename-readsname-windowsize.gz)')

opt_metric = make_option('-m', '--metric', type='choice',
        choices=['Mean counts', 'Frequency counts'],
        default='Frequency counts',
        help='Select the metric (note you will need to change your ylim accordingly)')

# sampling options
opt_cutoff = make_option('-k', '--cutoff', type='float', default = 0.05,
        help='Probability cutoff. Exclude genes if the probability of '\
        'the observed tag count is at most this value [default: %default]')
opt_top = make_option('-S', '--sample_top', type='int', default = None,
             help='Genes ranked to this number will be sampled (default is All)')
opt_stderr = make_option('-e', '--plot_stderr',
             action='store_true', help="Plot mean standard errors",
             default=False)

# optional plotting information
opts_legend1 = make_option('--legend1', help='Legend for sample 1')
opts_legend2 = make_option('--legend2', help='Legend for sample 2')

opt_fig_height = make_option('-H', '--fig_height', type='float',
        default=2.5*6, help='Figure height (cm) [default: %default]')

opt_fig_width = make_option('-W', '--fig_width', type='float', default=2.5*12,
   help='Figure width (cm) [default: %default]')
opt_bgcolor = make_option('-b', '--bgcolor', type='choice', default='black',
               help='Plot background color [default: %default]',
               choices=['black', 'white'])
opt_yrange = make_option('-y', '--ylim', default=None,
       help='comma separated minimum-maximum yaxis values (e.g. 0,3.5)')
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
opts_plot_filename = make_option('--plot_filename',
    default = None,
    help='Name of final plot file (must end with .pdf) [default: %default]')

# 
opt_test_run = make_option('-t', '--test_run',
             action='store_true', help="Test run, don't write output",
             default=False)

script_info['required_options'] = [opt_collection1, opts_legend1,
                    opt_collection2, opts_legend2, opt_metric, opt_yrange]

run_opts = [opt_test_run]
sampling_opts = [opt_cutoff, opt_top, opt_stderr]
save_opts = [opts_plot_filename]
plot_labels = [opts_title, opts_ylabel,
                opts_xlabel]
plot_dims = [opt_fig_height, opt_fig_width, opts_xgrid_locate,
            opts_ygrid_locate, opts_xlabel_interval, opts_ylabel_interval]
plot_colors = [opt_bgcolor, opts_vline_style, opts_vline_width,
        opts_xlabel_font, opts_ylabel_font]

script_info['optional_options'] = run_opts+sampling_opts+save_opts+\
        plot_labels+plot_dims+plot_colors

script_info['optional_options_groups'] = [('Run control', run_opts),
                                  ('Sampling', sampling_opts),
                                  ('Saving', save_opts),
                                  ('Plot text', plot_labels),
                                  ('Plot colours', plot_colors),
                                  ('Plot dimensions', plot_dims)
                                  ]


def plot_sample(plot, coll, calc_stat, x, title, xlabel, ylabel, color, label, stderr=False):
    if stderr:
        print '\tJackknifing the statistic and standard error'
        jk = JackknifeStats(coll.N, calc_stat)
        y = jk.JackknifedStat
        stderr = jk.StandardError
    else:
        stderr = None
        y = calc_stat(range(coll.N))
    plot(x, y=y, color=color, label=label, stderr=stderr, ylabel=ylabel, xlabel=xlabel, title=title)

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    
    if ',' not in opts.ylim:
        raise RuntimeError('ylim must be comma separated')
    
    ylim = map(float, opts.ylim.strip().split(','))
    
    print 'Loading counts data'
    data_collection1 = RegionCollection(filename=opts.collection1)
    window_size = data_collection1.info['args']['window_size']
    data_collection2 = RegionCollection(filename=opts.collection2)
    
    # filter both
    if opts.cutoff < 0 or opts.cutoff > 1:
        raise RuntimeError('The cutoff must be between 0 and 1')
    
    data_collection1 = data_collection1.filteredChebyshevUpper(opts.cutoff)
    data_collection2 = data_collection2.filteredChebyshevUpper(opts.cutoff)
    
    # make sure each collection consists ot the same genes
    shared_labels = set(data_collection1.labels) & \
                    set(data_collection2.labels)
    
    data_collection1 = data_collection1.filteredByLabel(shared_labels)
    data_collection2 = data_collection2.filteredByLabel(shared_labels)
    assert set(data_collection1.labels) == set(data_collection2.labels)
    
    if opts.sample_top is None:
        sample_top = data_collection1.N
    else:
        sample_top = opts.sample_top
    
    indices = range(sample_top)
    data_collection1 = data_collection1.take(indices)
    data_collection2 = data_collection2.take(indices)
    
    print 'Starting to plot'
    if opts.bgcolor == 'black':
        grid = {'color': 'w'}
        bgcolor = '0.1'
        vline_color = 'w'
    else:
        grid = {'color': 'k'}
        vline_color = 'k'
        bgcolor = '1.0'
    
    vline = dict(x=0, linewidth=opts.vline_width,
                   linestyle=opts.vline_style, color=vline_color)
    
    plot = PlottableSingle(height=opts.fig_height/2.5,
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
    
    if opts.metric == 'Mean counts':
        stat = averaged
    else:
        data_collection1 = data_collection1.asfreqs()
        data_collection2 = data_collection2.asfreqs()
        stat = summed
    
    plot_sample(plot, data_collection1, stat_maker(stat, data_collection1), x,
        opts.title, opts.xlabel, opts.ylabel, 'b', opts.legend1,
        opts.plot_stderr)
    plot_sample(plot, data_collection2, stat_maker(stat, data_collection2), x,
        opts.title, opts.xlabel, opts.ylabel, 'r', opts.legend2,
        opts.plot_stderr)
    
    plot.legend()
    plot.show()
    if opts.plot_filename and not opts.test_run:
        plot.savefig(opts.plot_filename)
    else:
        print opts.plot_filename
    

if __name__ == '__main__':
    main()

