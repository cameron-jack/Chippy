#!/usr/bin/env python
#coding: utf-8
from __future__ import division

import os, sys, glob

sys.path.extend(['..'])

import numpy
from math import sqrt
from chippy.util.command_args import Args
from chippy.core.collection import RegionCollection, column_sum, column_mean, stdev
from chippy.express.db_query import make_session, get_gene_ids
from chippy.draw.plottable import PlottableGroups
from chippy.draw.plot_data import PlotLine
from chippy.util.run_record import RunRecord

from chippy.util.util import create_path, dirname_or_default

__author__ = 'Gavin Huttley, Cameron Jack'
__copyright__ = 'Copyright 2011-2013, Gavin Huttley, Anuj Pahwa, Cameron Jack'
__credits__ = ['Gavin Huttley', 'Cameron Jack']
__license__ = 'GPL'
__maintainer__ = 'Cameron Jack'
__email__ = 'cameron.jack@anu.edu.au'
__status__ = 'Release'
__version__ = '0.2'

class Study(object):
    """ Specifies the RegionCollection associated with an expression
            data set. A div_study is marked as such.
    Members: collection (a RegionCollection), window_radius,
            collection_label
    Methods: filterByGenes, filterByCutoff, normaliseByBases,
            asPlotLines
    """
    def __init__(self, collection_fn, counts_func,
            *args, **kwargs):
        super(Study, self).__init__(*args, **kwargs)
        rr = RunRecord('Study')
        # Keep the source file name for labelling purposes
        fn = collection_fn.split('/')[-1].rstrip('.gz')
        self.collection_label = fn.replace('_', ' ')
        try:
            self.data_collection = RegionCollection(filename=collection_fn)
        except IOError:
            rr.dieOnCritical('Collection will not load', collection_fn)

        # Frequency normalized counts need to be converted
        if counts_func is column_sum:
            self.data_collection = self.data_collection.asfreqs()
        self.counts_func = counts_func

        # Get feature window radius
        try:
            self.window_radius =\
                    self.data_collection.info['args']['window_radius']
        except KeyError:
            self.window_radius = len(self.data_collection.counts[0])/2

    def filterByGenes(self, db_path, chrom=None, include_sample=None,
            exclude_sample = None):
        """ keep only results that match selected genes """
        rr = RunRecord('filterByGenes')
        if not include_sample and not exclude_sample:
            return

        rr.addInfo('Starting no. of genes', self.data_collection.N)

        session = make_session(db_path)
        if include_sample:
            include_sample = include_sample.split(':')[0].strip()
        if include_sample:
            exclude_sample = exclude_sample.split(':')[0].strip()

        filter_gene_ids = get_gene_ids(session, chrom=chrom,
                include_target=include_sample, exclude_target=exclude_sample)

        self.data_collection =\
                self.data_collection.filteredByLabel(filter_gene_ids)
        rr.addInfo('Remaining genes', self.data_collection.N)

        if self.data_collection is None or\
                self.data_collection.ranks.max() == 0:
            rr.dieOnCritical('No genes remaining after filtering', 'Failure' )

        # total_features used to normalise coloring
        total_features = self.data_collection.ranks.max()
        self.data_collection.ranks /= total_features

    def filterByCutoff(self, cutoff=None):
        """ keep only results that pass Chebyshev cutoff """
        rr = RunRecord('filterByCutoff')

        rr.addInfo('Starting no. of genes', self.data_collection.N)

        # exclude outlier genes using one-sided Chebyshev
        if cutoff is not None and cutoff != 0.0:
            try:
                cutoff = float(cutoff)
                if cutoff < 0.0 or cutoff >= 1.0:
                    rr.addError('Cutoff out of range', cutoff)
                    rr.addInfo('Cutoff set to default', 0.05)
                    cutoff = 0.05
            except ValueError:
                rr.addError('Cutoff not given as float', cutoff)
                rr.addInfo('Cutoff set to default', 0.05)
                cutoff = 0.05
            # Do Chebyshev filtering

            self.data_collection =\
                    self.data_collection.filteredChebyshevUpper(p=cutoff)
            rr.addInfo('Used Chebyshev filter cutoff', cutoff)
            rr.addInfo('No. genes after normalisation filter',
                    self.data_collection.N)
        else:
            rr.addInfo('Outlier cutoff filtering', 'Off')

        if self.data_collection is None or\
                self.data_collection.ranks.max() == 0:
            rr.dieOnCritical('No data after filtering', 'Failure')

        # total_features used to normalise coloring
        total_features = self.data_collection.ranks.max()
        self.data_collection.ranks /= total_features

    def normaliseByBases(self):
        # This requires 'base count' to be present in the collection
        rr = RunRecord('normaliseByBases')
        try:
            norm_bases = self.data_collection.info['args']['base count']
        except KeyError:
            rr.addError('Info field not found', 'base count')
            return

        rr.addInfo('normalising by RPMs', float(1000000/norm_bases))
        normalised_counts = []
        for c in self.data_collection.counts:
            c = c * 1000000 / norm_bases
            normalised_counts.append(c)
        self.data_collection.counts = normalised_counts

    def _groupAllGeneCounts(self):
        """ Group counts for all genes and return as a single PlotLine.
            Called by asPlotLines or _groupNGeneCounts().
            Returns a list.
        """
        rr = RunRecord('_groupAllGeneCounts')
        counts, ranks, se = self.data_collection.transformed(\
                counts_func=self.counts_func)
        if not len(counts):
            rr.dieOnCritical('No counts data in', 'Study.groupAllGeneCounts')

        # Always name single lines by their collection name
        label = self.collection_label
        plot_lines = [PlotLine(counts, ranks, label, study=label, stderr=se)]
        return plot_lines

    def _groupNoGeneCounts(self):
        """ Don't group counts. Simply return a PlotLine for each set of
            counts.
            Called by asPlotLines()
        """
        rr = RunRecord('_groupNoGeneCounts')
        counts = self.data_collection.counts
        ranks = self.data_collection.ranks
        labels = self.data_collection.labels
        plot_lines = []
        for c,r,l in zip(counts, ranks, labels):
            if self.counts_func == stdev:
                stdev_ = c.std()
                if stdev_ > 0:
                    c = (c - c.mean()) / stdev_
                    plot_lines.append(PlotLine(c, r , l,
                            study=self.collection_label))
            else:
                plot_lines.append(PlotLine(c, r , l,
                        study=self.collection_label))

        # If no data was returned default to groupAllCollectionCounts
        if not len(plot_lines):
            rr.dieOnCritical('No data in collection', 'Failure')

        # If a single line is created label it with the collection name
        if len(plot_lines) == 1:
            plot_lines[0].label = [self.collection_label]

        return plot_lines

    def _groupNGeneCounts(self, group_size):
        """ Group counts for N genes and return as PlotLines. Defaults to
            _groupAllGeneCounts() if group size is too large.
            Called by asPlotLines()
        """
        rr = RunRecord('_groupNGeneCounts')
        plot_lines = []
        for index, (c,r,l,se) in enumerate(self.data_collection.\
                iterTransformedGroups(group_size=group_size,
                counts_func=self.counts_func)):
            plot_lines.append(PlotLine(c, rank=index, label=l,
                    study=self.collection_label, stderr=se))

        # If no data was returned default to groupAllCollectionCounts
        if not len(plot_lines):
            rr.addWarning('Defaulting to ALL features. Not enough '+\
                    'features for group of size', group_size)
            plotLines = self.groupAllGeneCounts()
            return plotLines

        # If a single line is created label it with the collection name
        if len(plot_lines) == 1:
            plot_lines[0].label = [self.collection_label]

        return plot_lines

    def asPlotLines(self, studies, group_size, group_location):
        """ returns a list of PlotLine objects from this study """
        rr = RunRecord('asPlotLines')

        if type(group_size) is str and group_size.lower() == 'all':
            plot_lines= self._groupAllGeneCounts()
        elif type(group_size) is int:
            if group_size == 1:
                plot_lines = self._groupNoGeneCounts()
            else:
                plot_lines = self._groupNGeneCounts(group_size)
        else:
            rr.dieOnCritical('group_size, wrong type or value',
                    [type(group_size), group_size])

        if group_location:
            rr.addInfo('grouping genes from location', group_location)
            if group_location.lower() == 'top':
                plot_lines = [plot_lines[0]]
            elif group_location.lower() == 'middle':
                plot_lines = [plot_lines[len(plot_lines)/2]]
            elif group_location.lower() == 'bottom':
                plot_lines = [plot_lines[-1]]

        rr.addInfo('Plottable lines from study', len(plot_lines))
        return plot_lines

def interpret_group_size(str_group_size):
    """ group_size is provided as a string so special term 'All' may be
        given to refer to the grouping of all genes. """
    try:
        group_size = int(str_group_size)
    except ValueError:
        group_size = 'All'
    return group_size

def load_studies(collections, counts_func):
    """ load all collection data and apply filtering if needed """
    rr = RunRecord('load_studies')
    # Parse glob file names
    collection_files = collections
    dir_name = os.path.dirname(collection_files)
    base_name = os.path.basename(collection_files)
    collection_file_names = [os.path.join(dir_name, p)\
            for p in glob.glob1(dir_name, base_name)]
    collection_file_names.sort()

    window_radii = []
    studies = []
    # Load data from each file
    for collection_file in collection_file_names:
        study = Study(collection_file, counts_func)
        if study is None:
            rr.dieOnCritical('Could not load study', collection_file)
        else:
            studies.append(study)
            window_radii.append(study.window_radius)

    # Find max common windows size
    if not len(studies):
        rr.dieOnCritical('No valid data files', 'Failure')

    window_radius = int(min(window_radii))

    rr.addInfo('Max common window radius', window_radius)
    rr.addInfo('Total data collections', len(studies))

    return studies, window_radius

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

def set_counts_function(metric):
    """ Sets the feature counting metric function"""
    rr = RunRecord('set_counts_function')
    if metric.lower() == 'mean counts':
        counts_func = column_mean
        rr.addInfo('Counts metric set to', 'column_mean')
    elif metric.lower() == 'frequency counts':
        counts_func = column_sum
        rr.addInfo('Counts metric set to', 'column_sum')
    elif metric.lower() == 'standard deviation':
        counts_func = stdev
        rr.addInfo('Counts metric set to', 'stdev')
    else:
        rr.dieOnCritical('Invalid count metric', metric)
    return counts_func

def div_plots(plot_lines, div_study_name):
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
    if len(ranked_plot_lines) == 0:
        rr.dieOnCritical('No plot lines.', 'Same study as div plot?')

    out_lines = []
    for ranked_index, div_index in zip(ranked_plot_lines, dividing_plot_lines):
        if ranked_index != div_index:
            rr.dieOnCritical('Div and study plot lines do not match',
                    [len(ranked_plot_lines), len(dividing_plot_lines)])
        for line in ranked_plot_lines[ranked_index]:
            line.counts /= dividing_plot_lines[ranked_index].counts
            out_lines.append(line)

    return out_lines

def set_plot_colors(plot_lines, studies, bgcolor, grey_scale, cmap = None):
    rr = RunRecord('set_plot_colors')
    # Hack time: this is just for David's Cell-cycle plots
    if len(studies) == 3 and len(plot_lines) == 3:
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

    elif len(studies) > 1:
        # spectral, jet, hsv, gist_rainbow, rainbow - all decent options
        cmap = 'rainbow'

    elif len(studies) == 1 and grey_scale:
        # grey-scale spectrum. Since background is white or black we can't
        # have lines that are pure black or white.
        if bgcolor == 'black':
            for i, line in enumerate(plot_lines):
                col = (240/len(plot_lines)) + ((240/len(plot_lines))*i)
                r = col; g = col; b = col
                line.color = '#%02x%02x%02x' % (r, g, b)
        else: # white background
            for i, line in enumerate(plot_lines):
                col = 240 - ((240/len(plot_lines))*i)
                r = col; g = col; b = col
                line.color = '#%02x%02x%02x' % (r, g, b)

    elif len(studies) == 1:
        # coolwarm, RdBu, jet - all decent options
        cmap = 'RdBu'

    return plot_lines, cmap

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
req_args = ['collection', 'metric', 'plot_filename']
opt_args = ['ylim', 'fig_height', 'fig_width',
        'xgrid_lines', 'ygrid_lines', 'grid_off', 'xtick_interval',
        'ytick_interval', 'clean_plot', 'bgcolor', 'colorbar', 'title',
        'xlabel', 'ylabel', 'xfont_size', 'yfont_size', 'legend',
        'legend_font_size', 'vline_style', 'vline_width', 'grey_scale',
        'line_alpha', 'chrom', 'include_target', 'exclude_target', 'group_size',
        'group_location', 'top_features', 'smoothing', 'binning', 'cutoff',
        'plot_series', 'text_coords', 'test_run', 'version',
        'div', 'normalise_by_RPM', 'confidence_intervals']

script_info['args'] = Args(required_args=req_args, optional_args=opt_args,
    positional_args=pos_args)
script_info['required_options'] = script_info['args'].req_cogent_opts
script_info['optional_options'] = script_info['args'].opt_cogent_opts

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
    """
    rr = RunRecord('plot_centred_counts')
    rr.addCommands(sys.argv)
    args = script_info['args'].parse()

    # 1: Set feature counting metric
    counts_func = set_counts_function(args.metric)

    # 2: Load studies
    print 'Loading counts data'
    studies, window_radius = load_studies(args.collection, counts_func)

    # 3: Load divisor study if provided
    if args.div is not None:
        div_studies, div_window_radius = load_studies(args.div, counts_func)
        if div_window_radius == window_radius:
            print 'Windows match - using div study'
            studies.append(div_studies[0])
            div_name = div_studies[0].collection_label
        else:
            rr.dieOnCritical('Differing Data and Div window sizes',
                    [div_window_radius, window_radius])
    else:
        div_name = None

    # 4: Normalise counts if require
    if args.normalise_by_RPM:
        for study in studies:
            study.normaliseByBases()

    # 5: Specify genes of interest to direct study
    for study in studies:
        study.filterByGenes(args.db_path, chrom=args.chrom,
                include_sample=args.include_target,
                exclude_sample=args.exclude_target)

    # 6: Filter studies by statistical cutoff
    for study in studies:
        study.filterByCutoff(args.cutoff)

    # 7: Create plot lines for each study in studies
    group_size = interpret_group_size(args.group_size)
    plot_lines = []
    for study in studies:
        lines = study.asPlotLines(studies, group_size,
                args.group_location)
        for line in lines:
            plot_lines.append(line)

    rr.addInfo('Total number of lines from all studies', len(plot_lines))

    # 8: smooth and/or bin plot lines as required
    if args.binning and args.binning > 0:
        for line in plot_lines:
            line.applyBinning(args.binning)
        rr.addInfo('plot_centred_counts', 'lines binned to width',
                args.binning)

    if args.smoothing and args.smoothing > 0:
        for line in plot_lines:
            line.applySmoothing(args.smoothing)
        rr.addInfo('plot_centred_counts', 'lines smoothed to width',
                args.smoothing)

    # 9: Do plot division if required
    if div_name:
        plot_lines = div_plots(plot_lines, div_name)

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
            yaxis_lims=ylim, xaxis_lims=(-window_radius, window_radius),
            xy_tick_spaces=(args.xgrid_lines, args.ygrid_lines),
            xy_tick_intervals=(args.xtick_interval, args.ytick_interval),
            xy_label_fontsizes=(args.xfont_size, args.yfont_size),
            vline=vline, ioff=True, colorbar=args.colorbar,
            clean=args.clean_plot)
    
    x = numpy.arange(-window_radius, window_radius)

    # 11: set line colors
    cmap = None
    plot_lines, cmap = set_plot_colors(plot_lines, studies,\
            args.bgcolor, args.grey_scale, cmap=cmap)

    # 12: Create plot
    plot(x, y_series=None, plot_lines=plot_lines, color_series=None,
            series_labels=series_labels, filename_series=filename_series,
            label_coords=label_coords, cmap=cmap,
            alpha=args.line_alpha, xlabel=args.xlabel,
            ylabel=args.ylabel, title=args.title, colorbar=args.colorbar,
            labels=None, labels_size=args.legend_font_size,
            plot_CI=args.confidence_intervals)

    # 13: save plots
    # if series, create directory
    if args.plot_series and not args.test_run:
        set_up_series_plots_dir(args.plot_filename)

    if args.plot_filename and not args.test_run:
        if '.pdf' in args.plot_filename.lower():
            plot.savefig(args.plot_filename, image_format='pdf')
        else:
            plot.savefig(args.plot_filename+'.pdf', image_format='pdf')
    else:
        print args.plot_filename
    
    rr.display()
    plot.show()

if __name__ == '__main__':
    main()
