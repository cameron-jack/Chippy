#!/usr/bin/env python
#coding: utf-8
from __future__ import division

import os, sys, glob

sys.path.extend(['..'])

from matplotlib import cm, colors

import numpy
from chippy.util.command_args import Args
from chippy.core.collection import column_sum, column_mean, stdev
from chippy.draw.plottable import PlottableGroups
from chippy.util.run_record import RunRecord
from chippy.util.studies import RegionStudy

from chippy.util.util import create_path, dirname_or_default

__author__ = 'Cameron Jack, Gavin Huttley'
__copyright__ = 'Copyright 2011-2013, Gavin Huttley, Cameron Jack, Anuj Pahwa'
__credits__ = ['Cameron Jack', 'Gavin Huttley']
__license__ = 'GPL'
__maintainer__ = 'Cameron Jack'
__email__ = 'cameron.jack@anu.edu.au'
__status__ = 'Release'
__version__ = '0.2'

def load_studies(collections, counts_func):
    """
        Load all collection data as list, and apply filtering if needed.
        Return the studies plus their max common up- & down- stream
        window size.
    """
    rr = RunRecord('load_studies')

    # Parse glob file names
    if len(collections) == 0:
        rr.dieOnCritical('Number of provided collection files', 0)

    collection_fns = []
    for collection_file in collections:
        if '*' in collections[0]:
            dir_name = os.path.dirname(collection_file)
            base_name = os.path.basename(collection_file)
            glob_file_names = [os.path.join(dir_name, p)\
                    for p in glob.glob1(dir_name, base_name)]
            glob_file_names.sort()
            for f in glob_file_names:
                collection_fns.append(os.path.abspath(f))
        else:
            collection_fns.append(os.path.abspath(collection_file))

    windows_upstream = []
    windows_downstream = []
    studies = []
    # Load data from each file
    for collection_fn in collection_fns:
        study = RegionStudy(collection_fn, counts_func)
        if study is None:
            rr.dieOnCritical('Could not load study', collection_fn)
        else:
            studies.append(study)
            windows_upstream.append(study.window_upstream)
            windows_downstream.append(study.window_downstream)

    # Find max common windows size
    if not len(studies):
        rr.dieOnCritical('No valid data files', 'Failure')

    window_upstream = int(min(windows_upstream))
    window_downstream = int(min(windows_downstream))

    rr.addInfo('Max common upstream window size', window_upstream)
    rr.addInfo('Max common downstream window size', window_downstream)
    rr.addInfo('Total data collections', len(studies))

    return studies, window_upstream, window_downstream

def set_up_series_plots_dir(plot_filename):
    """ Create directory structure for series plots """
    rr = RunRecord('set_up_series_plot_dir')

    save_dir = dirname_or_default(plot_filename)
    basename = os.path.basename(plot_filename)

    plot_series_dir = os.path.join(save_dir,
        '%s-series' % basename[:basename.rfind('.')])
    create_path(plot_series_dir)
    rr.addInfo('Plotting as a series to', plot_series_dir)
    return plot_series_dir

def set_counts_function(counts_metric):
    """ Sets the feature counting metric function"""
    rr = RunRecord('set_counts_function')
    if counts_metric.lower() == 'mean':
        counts_func = column_mean
        rr.addInfo('Counts metric set to', 'column_mean')
    elif counts_metric.lower() == 'frequency':
        counts_func = column_sum
        rr.addInfo('Counts metric set to', 'column_sum')
    elif counts_metric.lower() == 'stdev':
        counts_func = stdev
        rr.addInfo('Counts metric set to', 'stdev')
    else:
        rr.dieOnCritical('Invalid count metric', counts_metric)
    return counts_func

def safe_line_division(sample_line, dividing_line):
    """ in case of divide-by-zero we need some robustness """
    rr = RunRecord('safe_line_division')
    try:
        sample_line /= dividing_line
    except ZeroDivisionError: # 0 counts at a base position
        min_count = 1
        for c in dividing_line:
            if 0 < c < min_count:
                min_count = c
        sample_line /= dividing_line
        rr.addWarning('Zero counts value seen. Setting zeros to',
                min_count)
    return sample_line

def div_plots(plot_lines, div_study_name, div_by='all'):
    """ Divides the counts values in plot_lines by those in the divisor
            lines """
    rr = RunRecord('div_plots')
    # build two matching sized dicts of lines indexed by rank
    ranked_plot_lines = {}
    dividing_plot_lines = {}
    for line in plot_lines:
        if line.study == div_study_name:
            dividing_plot_lines[line.rank] = line
        else:
            if not line.rank in ranked_plot_lines.keys():
                 ranked_plot_lines[line.rank] = []
            ranked_plot_lines[line.rank].append(line)

        # sanity check
        if len(ranked_plot_lines.keys()) == 0:
            rr.dieOnCritical('No plot lines.', 'Same study as div plot?')

    # default: divide scores line for line - maybe important if trends change
    out_lines = []
    if div_by.lower() == 'all':
        for ranked_index, div_index in zip(sorted(ranked_plot_lines.keys()),
                sorted(dividing_plot_lines.keys())):
            if ranked_index != div_index:
                rr.addCritical('Div and study plot lines do not match',
                        [len(ranked_plot_lines), len(dividing_plot_lines)])
                rr.addCritical('Differences in filtering may be responsible')
                rr.dieOnCritical('Try choosing a representative group instead',
                        ['mean-counts', 'median-counts',
                        'top-expr', 'median-expr', 'bottom-expr'])
            for line in ranked_plot_lines[ranked_index]:
                line.counts = safe_line_division(line.counts,
                        dividing_plot_lines[ranked_index].counts)
                out_lines.append(line)
    else: # divide all lines by the same counts array
        # build counts array to choose dividing counts from
        counts_lines = []
        div_counts = None
        for key in sorted(dividing_plot_lines.keys()):
            counts_lines.append(dividing_plot_lines[key].counts)
        counts_array = numpy.array(counts_lines)

        if div_by.lower() == 'mean-counts': # mean counts for all div
            div_counts = numpy.mean(counts_array, axis=0)
        elif div_by.lower() == 'median-counts': # median counts for all div
            div_counts = numpy.median(counts_array, axis=0)
        elif div_by.lower() == 'top-expr': # top expressed genes counts line from div
            div_counts = counts_array[0]
        elif div_by.lower() == 'median-expr': # median expressed genes counts line from div
            div_counts = counts_array[int(len(counts_array)/2)]
        elif div_by.lower() == 'bottom-expr': # bottom expressed genes counts line from div
            div_counts = counts_array[-1]
        else:
            rr.dieOnCritical('Unrecognised div_by type', div_by)

        for ranked_index in ranked_plot_lines:
            for line in ranked_plot_lines[ranked_index]:
                line.counts = safe_line_division(line.counts, div_counts)
                out_lines.append(line)

    return out_lines

def set_plot_colors(plot_lines, studies, div_name, bgcolor, grey_scale,
        restrict_colors=None):
    rr = RunRecord('set_plot_colors')
    min_col = 0
    max_col = 1
    if restrict_colors:
        col_limits = restrict_colors.split(',')
        min_col, max_col = (float(col_limits[0]), float(col_limits[1]))

    # Hack time: this is just for David's Cell-cycle plots
    num_studies = len(studies)
    if div_name: # one study doesn't count for coloring
        num_studies -= 1

    if num_studies == 3 and len(plot_lines) == 3:
        # We're going to use black/white, green and magenta, for G1, M, S
        if bgcolor == 'black':
            r = 255; g = 255; b = 255 # white
        else:
            r = 0; g = 0; b = 0 # black
        plot_lines[0].color = '#%02x%02x%02x' % (r, g, b)

        # green for M
        r = 0; g = 130; b = 0
        plot_lines[1].color = '#%02x%02x%02x' % (r, g, b)

        # Magenta for S
        r = 255; g = 0; b = 255
        plot_lines[2].color = '#%02x%02x%02x' % (r, g, b)

    elif num_studies > 1:
        # give each study a new colour and increase decrease alpha with rank

        # separate out plotlines per study
        per_study_lines = {}
        for line in plot_lines:
            study = line.study
            if study in per_study_lines.keys():
                per_study_lines[study].append(line)
            else:
                per_study_lines[study] = [line]

        plot_lines = []
        for i,s in zip(numpy.linspace(min_col, max_col, num_studies),
                per_study_lines.keys()):
            study_color = cm.jet(i) # cm.rainbow(i)
            rgba = colors.colorConverter.to_rgba(study_color)
            for study in studies:
                if s == study.collection_label:
                    for l,alpha in zip(sorted(per_study_lines[s],
                            key=lambda x: x.rank), numpy.linspace\
                            (0, 0.7, len(per_study_lines[s]))):
                        r, g, b, a = rgba
                        l.color = (r,g,b,a-alpha)
                        l.alpha = a-alpha
                        plot_lines.append(l)

    elif num_studies == 1 and grey_scale:
        # grey-scale spectrum. Since background is white or black we can't
        # have lines that are pure black or white.
        if bgcolor == 'black':
            for i, line in enumerate(plot_lines):
                # todo: replace this is with linspace and black/white
                col = (240/len(plot_lines)) + ((240/len(plot_lines))*i)
                r = col; g = col; b = col
                line.color = '#%02x%02x%02x' % (r, g, b)
        else: # white background
            for i, line in enumerate(plot_lines):
                col = 240 - ((240/len(plot_lines))*i)
                r = col; g = col; b = col
                line.color = '#%02x%02x%02x' % (r, g, b)

    elif num_studies == 1:
        # coolwarm, RdBu, jet - all decent options

        # make sure plots get coloured such that genes with low expr (large
        # number for rank get plotted blue
        plot_lines = sorted(plot_lines, key=lambda p: p.rank)

        # For our own plots use 0.05 and 0,85 for min_col, max_col
        for l,i in zip(plot_lines, numpy.linspace(min_col, max_col,
                len(plot_lines))):
            l.color = cm.RdBu(i)

    return plot_lines

script_info = {}
script_info['title'] = 'Plot read counts heat-mapped by gene expression'
script_info['script_description'] = 'Takes read counts that are centred on '\
        'on a gene feature such as TSS or intron-exon boundary, sorted '\
        'from high to low gene expression and makes a heat-mapped line plot.'
script_info['brief_description'] = 'Plots read counts around gene features'
script_info['version'] = __version__
script_info['authors'] = __author__
script_info['output_description']= 'Generates either a single pdf figure or '\
        'a series of pdfs that can be merged into a movie.'

pos_args = ['db_path']
req_args = ['collections', 'counts_metric']
opt_args = ['plot_filename', 'ylim', 'fig_height', 'fig_width',
        'xgrid_lines', 'ygrid_lines', 'grid_off', 'xtick_interval',
        'ytick_interval', 'offset_ticks', 'clean_plot', 'bgcolor', 'colorbar',
        'font', 'title', 'title_size', 'xlabel', 'ylabel', 'xfont_size',
        'yfont_size', 'legend',
        'legend_font_size', 'vline_style', 'vline_width', 'grey_scale',
        'line_width', 'restrict_colors', 'line_alpha', 'chrom',
        'include_targets', 'exclude_targets', 'group_size', 'group_location',
        'smoothing', 'binning', 'data_cutoff', 'line_cutoff', 'plot_series',
        'text_coords', 'test_run', 'version', 'div', 'div_by',
        'no_normalise', 'confidence_intervals',
        'write_genes_by_rank']

script_info['args'] = Args(required_args=req_args, optional_args=opt_args,
    positional_args=pos_args)
script_info['required_options'] = script_info['args'].getReqCogentOpts()
script_info['optional_options'] = script_info['args'].getOptCogentOpts()

def main():
    """ 1) Set counts_func
        2) Load studies
        3) Load divisor study if provided
        4) Normalise studies if required
        5) Set genes_of_interest
        6) Filter studies by genes_of_interest and statistical cutoff
        7) Create plotlines from studies
        8) Smooth or bin plotlines as required
        9) Do plot division (if required)
        10) Set basic plotting info
        11) Set lines colors as needed
        12) Create Plot
        13) Save Plot
        14) Save genes in plot to file in rank order (optional)
    """
    rr = RunRecord('plot_counts')
    rr.addCommands(sys.argv)
    args = script_info['args'].parse(window_title='Plot Counts')

    # 1: Set feature counting metric
    counts_func = set_counts_function(args.counts_metric)

    # 2: Load studies
    print 'Loading counts data'
    studies, window_upstream, window_downstream =\
            load_studies(args.collections, counts_func)

    # 3: Load divisor study if provided
    if args.div is not None:
        div_studies, div_window_upstream, div_window_downstream =\
                load_studies([args.div], counts_func)
        if div_window_upstream == window_upstream and \
                div_window_downstream == window_downstream:
            print 'Windows match - using div study'
            studies.append(div_studies[0])
            # alter name so we divide by the same study
            div_studies[0].collection_label += '_div'
            div_name = div_studies[0].collection_label
        else:
            rr.dieOnCritical('Differing Data and Div up/down-stream '+\
                    'window sizes',
                    [div_window_upstream, div_window_downstream,
                     window_upstream, window_downstream])
    else:
        div_name = None

    # 4: RPM Normalise counts by default
    if not args.no_normalise and args.counts_metric == 'mean':
        print 'Normalising by counts RPM'
        for study in studies:
            study.normaliseByRPM()

    # 5: Specify genes of interest to direct study
    for study in studies:
        study.filterByGenes(args.db_path, chrom=args.chrom,
                include_samples=args.include_targets,
                exclude_samples=args.exclude_targets)

    # 6: Filter all genes in studies by statistical cutoff
    if args.data_cutoff > 0.0:
        for study in studies:
            study.filterByCutoff(args.data_cutoff)

    # 7: Create plot lines for each study in studies
    try:
        group_size = int(args.group_size)
    except ValueError:
        group_size = 'All'

    plot_lines = []
    for study in studies:
        lines = study.asPlotLines(group_size, args.group_location,
                p=args.line_cutoff)

        for line in lines:
            plot_lines.append(line)

    rr.addInfo('Total number of lines from all studies', len(plot_lines))

    # 8: smooth and/or bin plot lines as required
    if args.binning and args.binning > 0:
        for line in plot_lines:
            line.applyBinning(args.binning)
        rr.addInfo('lines binned to width', args.binning)

    if args.smoothing and args.smoothing > 0:
        for line in plot_lines:
            line.applySmoothing(args.smoothing)
        rr.addInfo('lines smoothed to width', args.smoothing)

    # 9: Do plot division if required
    if div_name:
        plot_lines = div_plots(plot_lines, div_name, div_by=args.div_by)

    rr.addInfo('Total number of lines to plot', len(plot_lines))

    # 10: set basic plotting info
    ylim = None
    if args.ylim is not None:
        if ',' not in args.ylim:
            rr.dieOnCritical('ylim must be comma separated', ylim)
        ylim = map(float, args.ylim.strip().split(','))

    # if we have a plot series, create a directory to write plots
    if args.plot_series and not args.test_run:
        plot_series_dir = set_up_series_plots_dir(args.plot_filename)
        filename_series = []
    else:
        plot_series_dir = None
        filename_series = None
        series_labels = None
        label_coords = None
    
    print 'Prepping for plot'

    vline = dict(x=0, linewidth=args.vline_width,
            linestyle=args.vline_style, color='w')
    
    plot = PlottableGroups(height=args.fig_height/2.5,
            width=args.fig_width/2.5, bgcolor=args.bgcolor,
            grid_off=args.grid_off,
            yaxis_lims=ylim, xaxis_lims=(-window_upstream, window_downstream),
            xy_tick_spaces=(args.xgrid_lines, args.ygrid_lines),
            xy_tick_intervals=(args.xtick_interval, args.ytick_interval),
            offset_ticks=args.offset_ticks, linewidth=args.line_width,
            title_size=args.title_size, font=args.font,
            xy_label_fontsizes=(args.xfont_size, args.yfont_size),
            vline=vline, ioff=True, colorbar=args.colorbar,
            clean=args.clean_plot)
    
    x = numpy.arange(-window_upstream, window_downstream)

    # 11: set line colors
    plot_lines = set_plot_colors(plot_lines, studies,
            div_name, args.bgcolor, args.grey_scale,
            restrict_colors=args.restrict_colors)

    # 12: Create plot
    plot(x, plot_lines=plot_lines, filename_series=filename_series,
            xlabel=args.xlabel, ylabel=args.ylabel,
            title=args.title, colorbar=args.colorbar,
            labels_size=args.legend_font_size, show_legend=args.legend,
            plot_CI=args.confidence_intervals)

    # 13: Save plots
    # if series, create directory
    if args.plot_series and not args.test_run:
        set_up_series_plots_dir(args.plot_filename)

    if args.plot_filename and not args.test_run:
        plot_fn = args.plot_filename
        if '.pdf' in plot_fn.lower():
            plot.savefig(plot_fn, image_format='pdf')
        elif '.png' in plot_fn.lower():
            plot.savefig(plot_fn, image_format='png')
        elif '.jpg' in plot_fn.lower() or '.jpeg' in plot_fn.lower():
            plot.savefig(plot_fn, image_format='jpg')
        else:
            plot.savefig(plot_fn+'.pdf', image_format='pdf')
    else:
        plot.show()

    # 14: Save ENSEMBL gene ids by rank if requested
    if args.write_genes_by_rank:
        for study in studies:
            fn_parts = args.write_genes_by_rank.split('.')
            if len(fn_parts) > 1:
                f_ext = fn_parts[-1]
            else:
                f_ext = 'txt'
            fn = '.'.join(fn_parts[:-1]) + '_' +\
                    study.collection_label.replace(' ', '_') + '.' + f_ext
            with open(fn, 'w') as out:
                out.write('gene' + '\n')
                plot_lines.sort(key=lambda x: x.rank)
                for line in plot_lines:
                    for label in line.getLabelsAsList():
                        out.write(label + '\n')

    rr.display()

if __name__ == '__main__':
    main()
