from __future__ import division
from math import log10, floor, ceil

import os, sys, glob
sys.path.extend(['..', '../src'])

import numpy

from optparse import make_option
from cogent.util.misc import parse_command_line_parameters

from chippy.core.count_tags import centred_counts_for_genes,\
            centred_counts_external_genes
from chippy.core.collection import RegionCollection, column_sum, column_mean, column_stdev
from chippy.express import db_query
from chippy.draw.plottable import PlottableGroups
from chippy.ref.util import chroms
from chippy.util.run_record import RunRecord
from chippy.util.definition import LOG_DEBUG, LOG_INFO, LOG_WARNING, \
    LOG_ERROR, LOG_CRITICAL

from chippy.util.util import create_path, make_cl_command, just_filename, \
                    dirname_or_default

__author__ = 'Gavin Huttley, Cameron Jack'
__copyright__ = 'Copyright 2011, Gavin Huttley, Anuj Pahwa, Cameron Jack'
__credits__ = ['Gavin Huttley, Cameron Jack']
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

def _auto_grid_lines(ylim, test_run):
    """returns a float that is a 'round' looking number to use for the
            grid lines"""
    ymax = max(ylim)

    if ymax > 0:
        ypower = log10(ymax)

        if ypower < 0:
            rounding_places = 0 - int(floor(ypower))
            y_ceiling = float(ceil(ymax*(10**rounding_places))/(10**rounding_places))
            grid_lines = y_ceiling/10.0
        else:
            y_ceiling = ceil(ymax)
            if y_ceiling < 10:
                grid_lines = round(y_ceiling/10.0, 1)
            else:
                grid_lines = y_ceiling/10.0

    else:
        raise RuntimeError('Exiting: Maximum y-axis value meaningless: %e' % ymax)

    if test_run:
        print 'Setting Y-grid-line spacing: %e' % grid_lines
    return grid_lines

def _auto_yaxis(counts, ranks, test_run):
    """returns a list length 2 of y-axis limits, and matching grid_line value
    with 'round'-appearing numbers to make plots look pretty"""
    num_range = len(counts)
    ymaxs = []
    ymins = []
    for i in range(num_range):
        if ranks is not None:
            y = counts[i]
        else:
            y = counts

        ymaxs.append(max(y))
        ymins.append(min(y))

    ymax = max(ymaxs)
    ymin = min(ymins)

    ylim = (ymin, ymax)

    rounding_places = 1
    # For fractional counts then scale the rounding appropriately
    if ymax > 0:
        ypower = log10(ymax)
        if ypower < 0:
            rounding_places = 0 - int(floor(ypower))
            y_ceiling = float(ceil(ymax*(10**rounding_places))/(10**rounding_places))
            y_floor = float(floor(ymin*(10**rounding_places))/(10**rounding_places))
            grid_lines = y_ceiling/10.0
            ylim = (y_floor, y_ceiling)
        else:
            y_ceiling = ceil(ymax)
            y_floor = floor(ymin)
            if y_ceiling < 10:
                grid_lines = round(y_ceiling/10.0, 1)
            else:
                grid_lines = y_ceiling/10.0
                ylim = (y_floor,y_ceiling)
    elif ymax == 0:
        ylim = (0,1)
        grid_lines = 0.1
    else:
        raise RuntimeError('Exiting: Maximum y-axis value somehow negative: %e' % ymax)

    if test_run:
        print 'Y-max: %e, Y-min: %e' % (ymax, ymin)
        print 'Setting plot limits at Y-max: %e, Y-,min: %e' % (max(ylim), min(ylim))
        print 'Setting Y-grid-line spacing: %e' % grid_lines
    return ylim, grid_lines

def _filter_collection(data_collection, cutoff, external_sample, stable_ids, rr):
    # exclude outlier genes using one-sided Chebyshev
    if cutoff < 0 or cutoff > 1:
        raise RuntimeError('The cutoff must be between 0 and 1')

    rr.addMessage('plot_centred_counts._filter_collection', LOG_INFO,
        'Starting no. of genes', data_collection.N)
    if external_sample is None:
        data_collection = data_collection.filteredChebyshevUpper(p=cutoff)
        rr.addMessage('plot_centred_counts._filter_collection', LOG_INFO,
            'Used Chebyshev filter cutoff', cutoff)
        rr.addMessage('plot_centred_counts_filter_collection', LOG_INFO,
            'No. genes after normalisation filter', data_collection.N)

    if stable_ids is not None:
        data_collection = data_collection.filteredByLabel(stable_ids)
        rr.addMessage('plot_centred_counts', LOG_INFO,
            'Filtered by stable_ids', data_collection.N)

    total_gene = data_collection.ranks.max() # used to normalise colouring
    data_collection.ranks /= total_gene

    window_size = data_collection.info['args']['window_size']

    return data_collection, window_size, rr

def _group_genes(data_collection, group_size, labels, counts_func, topgenes, plot_series, rr):
    if group_size=='All':
        counts, ranks = data_collection.transformed(counts_func=counts_func)
        num_groups = 1
        counts = [counts]
        ranks = [ranks]
    else:
        counts = []
        ranks = []
        group_size = group_size
        group_size = int(group_size)
        for index, (c,r,l) in enumerate(data_collection.iterTransformedGroups(
                            group_size=group_size, counts_func=counts_func)):
            counts.append(c)
            ranks.append(r)
            if plot_series:
                labels.append('Group %d' % index)

        num_groups = len(counts)
        if num_groups == 0:
            counts, ranks = data_collection.transformed(counts_func=counts_func)
            num_groups = 1
            counts = [counts]
            ranks = [ranks]
            rr.addMessage('plot_centred_counts._group_genes', LOG_WARNING,
                'Defaulting to all genes. Not enough genes for group of size',
                group_size)

    if topgenes == True:
        num_groups = 1
        counts = list(counts.pop())
        ranks = None

    rr.addMessage('plot_centred_counts._group_genes', LOG_INFO,
        'Number of groups', num_groups)

    return counts, ranks, num_groups, labels, rr

if 'CHIPPY_DB' in os.environ:
    db_path = os.environ['CHIPPY_DB']
else:
    raise RuntimeError('You need to set an environment variable CHIPPY_DB '\
                       'that indicates where to find the database')

session = db_query.make_session('sqlite:///%s' % db_path)
samples = db_query.get_external_sample(session)

script_info = {}
script_info['title'] = 'Plot read counts heat-mapped by gene expression'
script_info['script_description'] = "Takes read counts that are centred on"\
    " on a gene TSS, sorted from high to low gene expression and makes a"\
    " heat-map plot."
script_info['version'] = __version__
script_info['authors'] = __author__
script_info['output_description']= "Generates either a single pdf figure or"\
    " a series of pdfs that can be merged into a movie."

# alternate option organisation

# essential source files
opt_collection = make_option('-s', '--collection',
  help='Path to the plottable data')

opt_metric = make_option('-m', '--metric', type='choice',
        choices=['Mean counts', 'Frequency counts', 'Standard deviation'],
        default='Frequency counts',
        help='Select the metric (note you will need to change your ylim '\
            'accordingly if providing via --ylim')
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

# optional sample choice information
opt_topgenes = make_option('--topgenes', action='store_true', default = False,
                    help='Plot only top genes ranked by expressed chromatin')

# optional plotting information
opt_cutoff = make_option('-k', '--cutoff', type='float', default = 0.05,
             help='Probability cutoff. Exclude genes if the probability of '\
             'the observed tag count is at most this value [default: %default]')
opt_fig_height = make_option('-H', '--fig_height', type='float', default=2.5*6,
   help='Figure height (cm) [default: %default]')

opt_fig_width = make_option('-W', '--fig_width', type='float', default=2.5*12,
   help='Figure width (cm) [default: %default]')
opt_bgcolor = make_option('-b', '--bgcolor', type='choice', default='black',
               help='Plot background color [default: %default]',
               choices=['black', 'white'])
opt_colorbar = make_option('--colorbar',
         action='store_true', help="Add colorbar to figure", default=False)
opt_yrange = make_option('-y', '--ylim', default=None,
       help='comma separated minimum-maximum yaxis values (e.g. 0,3.5)')
# Important note, grid_lines are an absolute scale!
opts_xgrid_locate = make_option('--xgrid_lines', type='float', default = 100,
                 help='major grid-line spacing on x-axis [default: %default]')
opts_ygrid_locate = make_option('--ygrid_lines', type='float', default = None,
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
    default = 'Normalized counts', help='Label for the y-axis [default: %default]')
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

script_info['required_options'] = [opt_collection, opt_metric]

run_opts = [opt_test_run]
sampling_opts = [opt_grp_size, opt_extern, opt_chroms, opt_cutoff, opt_topgenes]
save_opts = [opts_plot_filename]
series_opts = [opts_plotseries, opt_txt_coords]
plot_labels = [opts_title, opts_ylabel, opts_xlabel, opt_colorbar]
plot_dims = [opt_yrange, opt_fig_height, opt_fig_width, opts_xgrid_locate,
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
    
    rr = RunRecord()

    ylim = None
    if opts.ylim is not None:
        if ',' not in opts.ylim:
            raise RuntimeError('ylim must be comma separated')
        ylim = map(float, opts.ylim.strip().split(','))
    
    rr.addMessage('plot_centred_counts', LOG_INFO,
        'using metric', opts.metric)
    
    external_sample = get_sample_name(opts.external_sample)
    stable_ids = None
    if external_sample is not None:
        rr.addMessage('plot_centred_counts', LOG_INFO,
            'Using an external sample', external_sample)
        genes = db_query.get_external_genes(session, external_sample)
        stable_ids = [g.ensembl_id for g in genes]
    elif opts.chrom != 'All':
        rr.addMessage('plot_centred_counts', LOG_INFO,
            'Querying a single chromosome', opts.chrom)
        genes = db_query.get_genes(session, opts.chrom)
        stable_ids = [g.ensembl_id for g in genes]

    session.close()

    # if we have a plot series, we need to create a directory to dump the
    # files into
    if opts.plot_series and not opts.test_run:
        save_dir = dirname_or_default(opts.plot_filename)
        basename = os.path.basename(opts.plot_filename)
        # Should this be opts.plot_filename: ?
        #plot_filename = os.path.join(save_dir, basename)

        plot_series_dir = os.path.join(save_dir,
                        '%s-series' % basename[:basename.rfind('.')])
        create_path(plot_series_dir)
        rr.addMessage('plot_centred_counts', LOG_INFO,
            'Plotting as a series to', plot_series_dir)
        labels = []
        filename_series = []
    else:
        labels = None
        filename_series = None
        series_labels = None
        label_coords = None

    print 'Loading counts data'
    collection_files = opts.collection
    dir_name = os.path.dirname(collection_files)
    base_name = os.path.basename(collection_files)
    collection_file_names = [os.path.join(dir_name,
                p) for p in glob.glob1(dir_name, base_name)]

    window_size_set = []
    data_collection_set = []
    if opts.metric == 'Mean counts':
        print 'Calculating mean counts'
        counts_func = column_mean
        for collection_file in collection_file_names:
            data_collection = RegionCollection(filename=collection_file)
            # Filter genes for outliers and stableIDs
            data_collection, window_size, rr = _filter_collection(data_collection,
                    cutoff=opts.cutoff, external_sample=external_sample,
                    stable_ids=stable_ids, rr=rr)
            data_collection_set.append(data_collection)
            window_size_set.append(window_size)
            
    elif opts.metric == 'Frequency counts':
        print 'Calculating normalized frequency counts'
        counts_func = column_sum
        for collection_file in collection_file_names:
            data_collection = RegionCollection(filename=collection_file)
            data_collection = data_collection.asfreqs()
            # Filter genes for outliers and stableIDs
            data_collection, window_size, rr = _filter_collection(data_collection,
                    cutoff=opts.cutoff, external_sample=external_sample,
                    stable_ids=stable_ids, rr=rr)
            data_collection_set.append(data_collection)
            window_size_set.append(window_size)

    elif opts.metric == 'Standard deviation':
        print 'Calculating standard deviations of counts'
        counts_func = column_stdev
        for collection_file in collection_file_names:
            data_collection = RegionCollection(filename=collection_file)
            # Filter genes for outliers and stableIDs
            data_collection, window_size, rr = _filter_collection(data_collection,
                    cutoff=opts.cutoff, external_sample=external_sample,
                    stable_ids=stable_ids, rr=rr)
            data_collection_set.append(data_collection)
            window_size_set.append(window_size)

    else:
        print "--metric needs to be one of: 'Mean counts', 'Frequency counts', "\
              "or 'Standard deviation'"
        raise RuntimeError('Invalid metric choice')

    window_size = max(window_size_set)
    rr.addMessage('plot_centred_counts', LOG_INFO, 'Max window size', window_size)
    rr.addMessage('plot_centred_counts', LOG_INFO, 'Total data collections',
                  len(data_collection_set))

    # pool genes into groups
    count_set = []
    rank_set = []
    labels_set = []
    plottable_lines = 0 # total # of plotted lines
    for data_collection in data_collection_set:
        counts, ranks, num_groups, labels, rr = _group_genes(data_collection,
                group_size=opts.group_size, labels=labels,
                counts_func=counts_func, topgenes=opts.topgenes,
                plot_series=opts.plot_series, rr=rr)
        count_set.append(counts)
        rank_set.append(ranks)
        plottable_lines += num_groups
        labels_set.append(labels)

    rr.addMessage('plot_centred_counts', LOG_INFO,
        'Total number of plottable lines', plottable_lines)

    # reverse the counts and colour series so low color goes first
    if opts.plot_series:
        for labels in labels_set:
            label_coords = map(float, opts.text_coords.split(','))
            series_labels = list(reversed(labels))
            series_template = 'plot-%%.%sd.pdf' % len(str(len(counts)))
            filename_series = [os.path.join(plot_series_dir, series_template % i)
                            for i in range(len(series_labels))]
    
    print 'Prepping for plot'
    if opts.bgcolor == 'black':
        grid={'color': 'w'}
        bgcolor='0.1'
        vline_color='w'
    else:
        grid={'color': 'k'}
        vline_color='k'
        bgcolor='1.0'
    
    vline = dict(x=0, linewidth=opts.vline_width,
                   linestyle=opts.vline_style, color=vline_color)

    # auto-calculate y-min & y-max and/or y-tick-space, if required
    max_Ymax = None
    min_Ymin = None
    max_Ygrid_line = None
    if ylim is None:
        for counts, ranks in zip(count_set, rank_set):
            ylim, ygrid_line = _auto_yaxis(counts, ranks, opts.test_run)
            if (max_Ymax is None) or (max(ylim) > max_Ymax):
                max_Ymax = max(ylim)
                max_Ygrid_line = ygrid_line
            if (min_Ymin is None) or (min(ylim) < min_Ymin):
                min_Ymin = min(ylim)
        ylim = (min_Ymin, max_Ymax)
        opts.ygrid_lines = max_Ygrid_line
    else:
        if opts.ygrid_lines is None:
            opts.ygrid_lines = _auto_grid_lines(ylim, opts.test_run)

    maxY_str = '%e' % max(ylim)
    minY_str = '%e' % min(ylim)
    ygrid_line_str = '%e' % opts.ygrid_lines
    rr.addMessage('plot_centred_counts', LOG_INFO, 'Y-max plot limit',
                    maxY_str)
    rr.addMessage('plot_centred_counts', LOG_INFO, 'Y-min plot limit',
                    minY_str)
    rr.addMessage('plot_centred_counts', LOG_INFO, 'Y-grid-line spacing',
                    ygrid_line_str)
    
    plot = PlottableGroups(height=opts.fig_height/2.5,
        width=opts.fig_width/2.5,
        bgcolor=bgcolor, grid=grid,
        ylim=ylim, xlim=(-window_size, window_size),
        xtick_space=opts.xgrid_lines, ytick_space=opts.ygrid_lines,
        xtick_interval=opts.xlabel_interval,
        ytick_interval=opts.ylabel_interval,
        xlabel_fontsize=opts.xfontsize, ylabel_fontsize=opts.yfontsize,
        vline=vline, ioff=True, colorbar=opts.colorbar)
    
    x = numpy.arange(-window_size, window_size)

    all_ranks = []
    all_counts = []
    if len(data_collection_set) > 1:
        all_ranks = range(plottable_lines)
        for counts in count_set:
            for count in counts:
                all_counts.append(count)
        all_ranks = list(all_ranks)
        all_counts= list(reversed(all_counts))
    else:
        all_ranks = ranks
        all_counts = counts

    if opts.test_run:
        print 'Number of count sets: %d' % len(all_counts)
        if all_ranks is not None:
            print 'Number of rank sets: %d' % len(all_ranks)

    plot(x, y_series=all_counts, color_series=all_ranks, series_labels=series_labels,
        filename_series=filename_series, label_coords=label_coords,
        alpha=opts.line_alpha, xlabel=opts.xlabel,
        ylabel=opts.ylabel, title=opts.title, colorbar=opts.colorbar)
    
    if opts.plot_filename and not opts.test_run:
        plot.savefig(opts.plot_filename)
    else:
        print opts.plot_filename
    
    rr.display()
    plot.show()


if __name__ == '__main__':
    main()

