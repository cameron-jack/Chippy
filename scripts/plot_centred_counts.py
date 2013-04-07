from __future__ import division
from math import log10, floor, ceil

import os, sys, glob

sys.path.extend(['..'])

import numpy

from chippy.util.command_args import Args
from chippy.core.collection import RegionCollection, column_sum, column_mean, stdev
from chippy.express import db_query
from chippy.draw.plottable import PlottableGroups
from chippy.draw.util import smooth
from chippy.util.run_record import RunRecord
from chippy.util.definition import LOG_DEBUG, LOG_INFO, LOG_WARNING, \
    LOG_ERROR, LOG_CRITICAL

from chippy.util.util import create_path, dirname_or_default

__author__ = 'Gavin Huttley, Cameron Jack'
__copyright__ = 'Copyright 2011-2012, Gavin Huttley, Anuj Pahwa, Cameron Jack'
__credits__ = ['Gavin Huttley, Cameron Jack']
__license__ = 'GPL'
__maintainer__ = 'Cameron Jack'
__email__ = 'cameron.jack@anu.edu.au'
__status__ = 'alpha'
__version__ = '0.1'

def make_sample_choices(session):
    """returns the available choices for target gene samples"""
    samples = ['%s : %s' % (s.name, s.description)
        for s in db_query.get_target_samples(session)]
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
        if type(counts[i]) == numpy.float64:
            ymaxs.append(counts[i])
            ymins.append(counts[i])
        else:
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

def _filter_collection(data_collection, cutoff, target_sample, stable_ids, rr):
    # exclude outlier genes using one-sided Chebyshev
    if cutoff < 0 or cutoff > 1:
        raise RuntimeError('The cutoff must be between 0 and 1')

    rr.addMessage('plot_centred_counts._filter_collection', LOG_INFO,
        'Starting no. of genes', data_collection.N)
    if target_sample is None:
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

    try:
        window_size = data_collection.info['args']['window_size']        
    except KeyError:
        print 'No info tags'
        window_size = len(data_collection.counts[0])/2  
    
    return data_collection, window_size, rr

def _group_genes(data_collection, group_size, labels, counts_func, top_features, plot_series, rr):

    if group_size=='All':
        counts, ranks = data_collection.transformed(counts_func=counts_func)
        num_groups = 1
        counts = [counts]
        ranks = [ranks]
        labels_set = [labels]
    else:
        counts = []
        ranks = []
        labels_set = []
        group_size = group_size
        group_size = int(group_size)
        for index, (c,r,l) in enumerate(data_collection.iterTransformedGroups(
                            group_size=group_size, counts_func=counts_func)):
            counts.append(c)
            ranks.append(r)
            labels_set.append(labels)

            if top_features:
                num_groups = 1
                return counts, ranks, num_groups, labels_set, rr

            if plot_series:
                labels_set.append('Group %d' % index)

        num_groups = len(counts)
        if num_groups == 0: # default to 1 group
            counts, ranks = data_collection.transformed(counts_func=counts_func)
            num_groups = 1
            counts = [counts]
            ranks = [ranks]
            labels_set = [labels]
            rr.addMessage('plot_centred_counts._group_genes', LOG_WARNING,
                'Defaulting to all genes. Not enough genes for group of size',
                group_size)

    rr.addMessage('plot_centred_counts._group_genes', LOG_INFO,
        'Number of groups', num_groups)

    return counts, ranks, num_groups, labels_set, rr

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
        'legend_font_size', 'vline_style', 'vline_width',
        'line_alpha', 'chrom', 'external_sample', 'group_size',
        'group_location', 'top_features', 'smoothing', 'binning', 'cutoff',
        'plot_series', 'text_coords', 'test_run', 'version',
        'div', 'normalise_tags1', 'normalise_tags2', 'normalise_tags3',
        'normalise_by_RPM']

script_info['args'] = Args(required_args=req_args, optional_args=opt_args,
    positional_args=pos_args)
script_info['required_options'] = script_info['args'].req_cogent_opts
script_info['optional_options'] = script_info['args'].opt_cogent_opts

def main():
    rr = RunRecord()
    args = script_info['args'].parse()

    ylim = None
    if args.ylim is not None:
        if ',' not in args.ylim:
            raise RuntimeError('ylim must be comma separated')
        ylim = map(float, args.ylim.strip().split(','))
    
    rr.addMessage('plot_centred_counts', LOG_INFO,
        'using metric', args.metric)

    target_sample=None
    #target_sample = get_sample_name(args.target_sample)
    stable_ids = None
    #if target_sample is not None:
    #    rr.addMessage('plot_centred_counts', LOG_INFO,
    #        'Using an target sample', target_sample)
    #    genes = db_query.get_target_genes(session, target_sample)
    #    stable_ids = [g.ensembl_id for g in genes]
    #elif args.chrom != 'All':
    #    rr.addMessage('plot_centred_counts', LOG_INFO,
    #        'Querying a single chromosome', args.chrom)
    #    genes = db_query.get_genes(session, args.chrom)
    #    stable_ids = [g.ensembl_id for g in genes]

    # if we have a plot series, we need to create a directory to dump the
    # files into
    if args.plot_series and not args.test_run:
        save_dir = dirname_or_default(args.plot_filename)
        basename = os.path.basename(args.plot_filename)

        plot_series_dir = os.path.join(save_dir,
                        '%s-series' % basename[:basename.rfind('.')])
        create_path(plot_series_dir)
        rr.addMessage('plot_centred_counts', LOG_INFO,
            'Plotting as a series to', plot_series_dir)
        filename_series = []
    else:
        filename_series = None
        series_labels = None
        label_coords = None


    print 'Loading counts data'
    collection_files = args.collection
    dir_name = os.path.dirname(collection_files)
    base_name = os.path.basename(collection_files)
    collection_file_names = [os.path.join(dir_name,
                p) for p in glob.glob1(dir_name, base_name)]
    collection_file_names.sort()
    filenames_set = []
    for file in collection_file_names:
        file_parts = file.split('/')
        file = file_parts[-1]
        file = file.rstrip('.gz')
        file = file.replace('_', ' ')
        filenames_set.append(file)

    print 'collection filename set size: %d' % len(collection_file_names)

    window_size_set = []
    data_collection_set = []
    if args.metric == 'Mean counts':
        print 'Collating mean counts'
        counts_func = column_mean
        for collection_file in collection_file_names:
            data_collection = RegionCollection(filename=collection_file)
            # Filter genes for outliers and stableIDs
            data_collection, window_size, rr = _filter_collection(data_collection,
                    cutoff=args.cutoff, target_sample=target_sample,
                    stable_ids=stable_ids, rr=rr)
            data_collection_set.append(data_collection)
            window_size_set.append(window_size)

    elif args.metric == 'Frequency counts':
        print 'Collating normalized frequency counts'
        counts_func = column_sum
        for collection_file in collection_file_names:
            data_collection = RegionCollection(filename=collection_file)
            data_collection = data_collection.asfreqs()
            # Filter genes for outliers and stableIDs
            data_collection, window_size, rr = _filter_collection(data_collection,
                    cutoff=args.cutoff, target_sample=target_sample,
                    stable_ids=stable_ids, rr=rr)
            data_collection_set.append(data_collection)
            window_size_set.append(window_size)

    elif args.metric == 'Standard deviation':
        print 'Collating standard deviations of counts'
        counts_func = stdev
        for collection_file in collection_file_names:
            data_collection = RegionCollection(filename=collection_file)
            # Filter genes for outliers and stableIDs
            data_collection, window_size, rr = _filter_collection(data_collection,
                    cutoff=args.cutoff, target_sample=target_sample,
                    stable_ids=stable_ids, rr=rr)
            data_collection_set.append(data_collection)
            window_size_set.append(window_size)

    else:
        print "--metric needs to be one of: 'Mean counts', 'Frequency counts', "\
              "or 'Standard deviation'"
        raise RuntimeError('Invalid metric choice')

    if len(window_size_set) == 0:
        raise RuntimeError('No valid data files loaded')

    window_size = min(window_size_set)
    rr.addMessage('plot_centred_counts', LOG_INFO, 'Max window size', window_size)
    rr.addMessage('plot_centred_counts', LOG_INFO, 'Total data collections',
                  len(data_collection_set))

    if args.group_size.lower() == 'all':
        group_size = 'All'
    else:
        try:
            group_size = int(args.group_size)
        except ValueError:
            print ('Invalid group size: ' + args.group_size + '. Defaulting to all genes.\n')
            group_size = 'All'

    # pool genes into groups
    count_set = []
    rank_set = []
    labels_set = []
    plottable_lines = 0 # total # of plotted lines
    iteration = 0
    for dc_index, data_collection in enumerate(data_collection_set):
        counts, ranks, num_groups, labels, rr = _group_genes(data_collection,
                group_size=group_size, labels=filenames_set[iteration],
                counts_func=counts_func, top_features=args.top_features,
                plot_series=args.plot_series, rr=rr)

        if args.smoothing > 0:
            smoothed_counts = []
            for c in counts:
                c = smooth(c, args.smoothing)              
                smoothed_counts.append(c)
            counts = smoothed_counts

        if group_size == 'All':
            genes_per_group = data_collection.N
        else:
            genes_per_group = group_size

        # Calculate normalised tags per million mapped reads (RPM)
        # Since we have base counts instead of tag counts we should normalise
        # by total base counts
        if args.metric.lower() == 'mean counts':
            if args.normalise_by_RPM:
                normalised_counts = []
                norm_bases = data_collection.info['args']['base count']
                for c in counts:
                    c = c * 1000000 / norm_bases
                    normalised_counts.append(c)
                counts = normalised_counts
            else:
                if args.normalise_tags1 is not None and dc_index == 0:
                    normalised_counts = []
                    norm_tags = args.normalise_tags1
                    for c in counts:
                        # Which is better, per line or per gene normalisation?
                        #c = c * genes_per_group * 1000000 / norm_tags
                        # ^- this is per group/line normalisation
                        c = c * 1000000 / norm_tags # This is per gene normalisation
                        normalised_counts.append(c)
                    counts = normalised_counts

                if args.normalise_tags2 is not None and dc_index == 1:
                    # Calculate normalised per million mapped reads (RPM)
                    # Hack for 3 data sets, G1, M, S in that order
                    normalised_counts = []
                    norm_tags = args.normalise_tags2
                    for c in counts:
                        #c = c * genes_per_group * 1000000 / norm_tags
                        c = c * 1000000 / norm_tags
                        normalised_counts.append(c)
                    counts = normalised_counts

                if args.normalise_tags3 is not None and dc_index == 2:
                    # Calculate normalised per million mapped reads (RPM)
                    # Hack for 3 data sets, G1, M, S in that order
                    normalised_counts = []
                    norm_tags = args.normalise_tags3
                    for c in counts:
                        #c = c * genes_per_group * 1000000 / norm_tags
                        c = c * 1000000 / norm_tags
                        normalised_counts.append(c)
                    counts = normalised_counts

        count_set.append(counts)
        rank_set.append(ranks)
        plottable_lines += num_groups
        for label in labels:
            labels_set.append(label)
        iteration += 1

    if args.div and args.top_features and len(count_set) == 2:
        # divide one set of counts by the other.
        # top100 genes is essential to keep plots comparable
        div_c = []
        if args.div == 1:
            # the first set is the divisor
            denominator_counts = count_set[0]
            numerator_counts = count_set[1]
        else:
            # the second set is the divisor
            denominator_counts = count_set[1]
            numerator_counts = count_set[0]

        d_counts = denominator_counts[0]
        n_counts = numerator_counts[0]

        for i in range(d_counts.size):
            try:
                div_c.append((float(n_counts[i])/float(d_counts[i])))
            except ZeroDivisionError:
                div_c.append(1.0)

        div_counts = numpy.array(div_c)
        counts = [div_counts]
        count_set = [counts]
        rank_set = [rank_set[0]]
        labels_set = [labels_set[0]]


    rr.addMessage('plot_centred_counts', LOG_INFO,
        'Total number of plottable lines', plottable_lines)

    # reverse the counts and colour series so low color goes first
    if args.plot_series:
        for labels in labels_set:
            label_coords = map(float, args.text_coords.split(','))
            series_labels = list(reversed(labels))
            series_template = 'plot-%%.%sd.pdf' % len(str(len(counts)))
            filename_series = [os.path.join(plot_series_dir, series_template % i)
                            for i in range(len(series_labels))]
    
    print 'Prepping for plot'
    if args.bgcolor == 'black':
        if args.grid_off is True:
            grid=False
            vline_color='k'
        else:
            grid={'color': 'w'}
            vline_color='w'
        bgcolor='0.1'
    else:
        if args.grid_off is True:
            grid=False
            vline_color='w'
        else:
            grid={'color': 'k'}
            vline_color='k'
        bgcolor='1.0'
    
    vline = dict(x=0, linewidth=args.vline_width,
                   linestyle=args.vline_style, color=vline_color)

    # auto-calculate y-min & y-max and/or y-tick-space, if required
    max_Ymax = None
    min_Ymin = None
    max_Ygrid_line = None
    if ylim is None:
        for counts, ranks in zip(count_set, rank_set):
            ylim, ygrid_line = _auto_yaxis(counts, ranks, args.test_run)
            if (max_Ymax is None) or (max(ylim) > max_Ymax):
                max_Ymax = max(ylim)
                max_Ygrid_line = ygrid_line
            if (min_Ymin is None) or (min(ylim) < min_Ymin):
                min_Ymin = min(ylim)
        ylim = (min_Ymin, max_Ymax)
        args.ygrid_lines = max_Ygrid_line
    else:
        if args.ygrid_lines is None:
            args.ygrid_lines = _auto_grid_lines(ylim, args.test_run)

    maxY_str = '%e' % max(ylim)
    minY_str = '%e' % min(ylim)
    ygrid_line_str = '%e' % args.ygrid_lines
    rr.addMessage('plot_centred_counts', LOG_INFO, 'Y-max plot limit',
            maxY_str)
    rr.addMessage('plot_centred_counts', LOG_INFO, 'Y-min plot limit',
            minY_str)
    rr.addMessage('plot_centred_counts', LOG_INFO, 'Y-grid-line spacing',
            ygrid_line_str)

    # Rather than have everything that follows simply dump into
    # PlottableGroups, it might be better to have multiple calls
    # to PlottableSingle
    
    plot = PlottableGroups(height=args.fig_height/2.5,
            width=args.fig_width/2.5, bgcolor=bgcolor, grid=grid,
            ylim=ylim, xlim=(-window_size, window_size),
            xtick_space=args.xgrid_lines, ytick_space=args.ygrid_lines,
            xtick_interval=args.xtick_interval,
            ytick_interval=args.ytick_interval,
            xlabel_fontsize=args.xfont_size, ylabel_fontsize=args.yfont_size,
            vline=vline, ioff=True, colorbar=args.colorbar,
            clean=args.clean_plot)
    
    x = numpy.arange(-window_size, window_size)

    all_ranks = []
    all_counts = []
    if len(count_set) > 1:
        all_ranks = range(plottable_lines)
        for counts in count_set:
            for count in counts:
                counts = list(reversed(counts))
                all_counts.append(count)
        all_ranks = list(all_ranks)
        all_counts= list(all_counts)
    else:
        counts = list(reversed(counts))
        ranks = list(reversed(ranks))

    if args.test_run:
        print 'Number of count sets: %d' % len(all_counts)
        if all_ranks is not None:
            print 'Number of rank sets: %d' % len(all_ranks)

    # Hack time: this is just for David's Cell-cycle plots
    if len(count_set) == 3:
        colour_range = []
        # We're going to use black/white, green and magenta, for G1, M, S
        if args.bgcolor == 'black':
            r = 255
            g = 255
            b = 255
        else:
            r = 0
            g = 0
            b = 0
        colour = '#%02x%02x%02x' % (r, g, b)
        colour_range.append((colour))
        # green for M
        r = 0
        g = 130
        b = 0
        colour = '#%02x%02x%02x' % (r, g, b)
        colour_range.append((colour))
        # Magenta for S
        r = 255
        g = 0
        b = 255
        colour = '#%02x%02x%02x' % (r, g, b)
        colour_range.append((colour))

        if not args.legend:
            labels_set = None

        plot(x, y_series=all_counts, color_series=colour_range,
                series_labels=series_labels, filename_series=filename_series,
                label_coords=label_coords, alpha=args.line_alpha,
                xlabel=args.xlabel, ylabel=args.ylabel, title=args.title,
                colorbar=args.colorbar, labels=labels_set,
                labels_size=args.legend_font_size)

    elif len(count_set) > 1:
        # spread colours almost evenly throughout the 256^3 colour-space
        colour_range = []
        halfway = floor(len(all_counts)/2)
        for i in range(len(all_counts)):
            if i < halfway:
                b = int(255 - ((255/halfway)*i))
                g = int(i * 255 / halfway)
                r = 0

            elif i == halfway:
                r = 0
                g = 255
                b = 0

            else:
                b = 0
                g = 255 - int((255/(halfway)) * (i-halfway))
                r = int((i-halfway)*(255/halfway))

            colour = '#%02x%02x%02x' % (r, g, b)
            if args.test_run:
                print 'Count colour for set %d is:' % i
                print colour

            colour_range.append(colour)

        if not args.legend:
            labels_set = None
        plot(x, y_series=all_counts, color_series=colour_range, series_labels=series_labels,
            filename_series=filename_series, label_coords=label_coords,
            alpha=args.line_alpha, xlabel=args.xlabel,
            ylabel=args.ylabel, title=args.title, colorbar=args.colorbar,
            labels=labels_set, labels_size=args.legend_size)
    else:
        if not args.legend:
            labels_set = None

        if len(counts) == 4: #quartile plot
            # do a b&w colour scheme
            color_range = []
            if args.bgcolor == 'black':
                for i in xrange(4):
                    r = 60 + (60*i)
                    g = 60 + (60*i)
                    b = 60 + (60*i)
                    color = '#%02x%02x%02x' % (r, g, b)
                    color_range.append(color)
            else: # white background
                for i in xrange(4):
                    r = 180 - (60*i)
                    g = 180 - (60*i)
                    b = 180 - (60*i)
                    color = '#%02x%02x%02x' % (r, g, b)
                    color_range.append(color)

            plot(x, y_series=counts, color_series=color_range, series_labels=series_labels,
                filename_series=filename_series, label_coords=label_coords,
                alpha=args.line_alpha, xlabel=args.xlabel,
                ylabel=args.ylabel, title=args.title, colorbar=args.colorbar,
                labels=labels_set, labels_size=args.legend_font_size)

        else:
            plot(x, y_series=counts, color_series=ranks, series_labels=series_labels,
                filename_series=filename_series, label_coords=label_coords,
                alpha=args.line_alpha, xlabel=args.xlabel,
                ylabel=args.ylabel, title=args.title, colorbar=args.colorbar,
                labels=labels_set, labels_size=args.legend_font_size)
    
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

