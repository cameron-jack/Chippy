from __future__ import division
from math import log10, floor, ceil

import sys
sys.path.extend(['..', '../src'])

from chippy.util.run_record import RunRecord
from matplotlib import pyplot, rcParams
from chippy.util.command_args import Args
from chippy.util.studies import MatchedStudy

__author__ = 'Cameron Jack'
__copyright__ = 'Copyright 2011-2013, Gavin Huttley, Cameron Jack, Anuj Pahwa'
__credits__ = ['Cameron Jack']
__license__ = 'GPL'
__maintainer__ = 'Cameron Jack'
__email__ = 'cameron.jack@anu.edu.au'
__status__ = 'pre-release'
__version__ = '0.2'

# Describe the application
script_info = {}
script_info['title'] = 'Plot gene mapped chromatin counts or ranks vs' \
        ' expression by rank or score'
script_info['script_description'] = 'Comparative dot or line plots of ' \
        'expression or mapped reads against the same or different data, ' \
        'ranked or unranked, for any number of studies. The number of '\
        'samples must be the same for each axis. '
script_info['usage'] = 'You can exclude genes '\
        'or force their inclusion by first using gene_overlap.py to ' \
        'generate a gene list and upload it to the DB with ' \
        'add_expression.py, then use --iX or --eX to Include or Exclude '\
        'from the corresponding --sX Sample.'
script_info['version'] = __version__
script_info['authors'] = __author__
script_info['output_description']= 'PNG/PDF histogram, line or dot plot'
script_info['help_on_no_arguments'] = True

pos_args = ['db_path']
req_args = ['sample', 'collection']
opt_args = ['x_axis_type', 'group_location', 'include_target',
            'fig_height', 'fig_width', 'region_feature',
            'exclude_target', 'plot_filename', 'counts_is_ranks',
            'expr_is_ranks', 'x_axis_is_log', 'y_axis_is_log']

script_info['args'] = Args(required_args=req_args, optional_args=opt_args,
    positional_args=pos_args)
script_info['required_options'] = script_info['args'].req_cogent_opts
script_info['optional_options'] = script_info['args'].opt_cogent_opts

class FigureDetails:
    """ This should be replaced with Plottable, which should also have a
    set_plot_options function so that options for setting title and such
    are standard and don't need to be reimplemented """

    def __init__(self, x_size=5, y_size=3, title=None, x_text=None,
            y_text=None):
        self.x_size = x_size
        self.y_size = y_size
        self.title = title
        self.x_text = x_text
        self.y_text = y_text

def make_plot(plot_points, plot_type='dot', plot_fn=None,
        output_type=None, fig_details=None, x_axis_is_log=False,
        y_axis_is_log=False, x_axis_type='expression', expr_is_ranks=False,
        counts_is_ranks=False):
    """
        Display plot_points as either dot or line plot.
        The points are always sorted by x-axis.
    """
    rr = RunRecord('make_plot')

    if not fig_details:
        fig_details = FigureDetails()

    fig = pyplot.figure(figsize=(fig_details.x_size, fig_details.y_size))
    ax = fig.add_subplot(111)
    pyplot.title(fig_details.title)
    pyplot.ylabel(fig_details.y_text)
    pyplot.xlabel(fig_details.x_text)
    rcParams['xtick.direction'] = 'out'
    rcParams['ytick.direction'] = 'out'
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    rr.addInfo('Output plot', plot_type)

    x_vals = []
    y_vals = []

    plot_points.sort()

    for point in plot_points:
    # transform plot points if required
        if x_axis_is_log:
            x_vals.append(point.get_logX())
        else:
            x_vals.append(point.x)

        if y_axis_is_log:
            y_vals.append(point.get_logY())
        else:
            y_vals.append(point.y)

    if plot_type == 'dot':
        ax.plot(x_vals, y_vals, 'o')
    elif plot_type == 'line':
        ax.plot(x_vals, y_vals)
    else:
       rr.addInfo('plot type not supported', plot_type)

    # If either axis is ranks (not scores) then we need to reverse the order
    # of the axes so we go from low to high
    if (x_axis_type.lower() == 'expression' and expr_is_ranks) or\
            (x_axis_type.lower() == 'counts' and counts_is_ranks):
        pyplot.gca().invert_xaxis()

    if (x_axis_type.lower() == 'expression' and counts_is_ranks) or\
       (x_axis_type.lower() == 'counts' and expr_is_ranks):
        pyplot.gca().invert_yaxis()

    if plot_fn:
        if plot_fn.lower().endswith('.pdf'):
            pyplot.savefig(plot_fn, format='pdf')
        elif plot_fn.lower().endswith('.jpg') or \
                plot_fn.lower().endswith('.jpeg'):
            pyplot.savefig(plot_fn, format='jpg')
        elif plot_fn.lower().endswith('.png'):
            pyplot.savefig(plot_fn, format='png')
        else:
            rr.addWarning('Unrecognised plot file extension. Using', 'PNG')
            pyplot.savefig(plot_fn + '.png', format='png')
    else:
        pyplot.show()

def main():
    """
        Comparative plots of count or rank data for chromatin or expression.
    """

    rr = RunRecord('counts_vs_expr')
    rr.addCommands(sys.argv)
    args = script_info['args'].parse(use_scrollbars=True,
        use_save_load_button=True,
        window_title='Counts vs Expression Plots')

    # Load all required data
    print 'Loading expression and counts data'
    matched_studies = MatchedStudy(args.sample, args.collection,
            args.db_path, args.region_feature,
            include_target=args.include_target,
            exclude_target=args.exclude_target)

    print 'Creating plot points'
    plot_points = matched_studies.get_matched_genes_as_xy_plotpoints(\
            args.x_axis_type, args.expr_is_ranks, args.counts_is_ranks)

    #fig = FigureDetails(x_size=args.fig_width, y_size=args.fig_height,
    #        title=args.sample + ' vs ' + args.collection)

    fig = FigureDetails(x_size=8, y_size=6,
        title=args.sample + ' vs ' + args.collection)

    if args.x_axis_type.lower() == 'expression':
        fig.x_text = 'Expression'
        fig.y_text = 'Counts'
        if args.counts_is_ranks:
            fig.y_text += ' Ranks'
        if args.expr_is_ranks:
            fig.x_text += ' Ranks'
    else:
        fig.y_text = 'Expression'
        fig.x_text = 'Counts'
        if args.counts_is_ranks:
            fig.x_text += ' Ranks'
        if args.expr_is_ranks:
            fig.y_text += ' Ranks'

    if args.x_axis_is_log:
        fig.x_text += ' (log base 2)'
    if args.y_axis_is_log:
        fig.y_text += ' (log base 2)'

    make_plot(plot_points, plot_fn=args.plot_filename,
            fig_details=fig, x_axis_is_log=args.x_axis_is_log,
            y_axis_is_log=args.y_axis_is_log, x_axis_type=args.x_axis_type,
            counts_is_ranks=args.counts_is_ranks,
            expr_is_ranks=args.expr_is_ranks)

    rr.display()

if __name__ == '__main__':
    main()

