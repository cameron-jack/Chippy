from __future__ import division

import sys
sys.path.extend(['..'])

import numpy

from chippy.util.run_record import RunRecord
from matplotlib import pyplot, rcParams
from chippy.draw.plottable import FigureDetails
from chippy.util.command_args import Args
from chippy.util.studies import CountsStudy

__author__ = 'Cameron Jack'
__copyright__ = 'Copyright 2011-2014, Gavin Huttley, Anuj Pahwa, Cameron Jack'
__credits__ = ['Cameron Jack']
__license__ = 'GPL'
__maintainer__ = 'Cameron Jack'
__email__ = 'cameron.jack@anu.edu.au'
__status__ = 'pre-release'
__version__ = '1.0'

script_info = {}
script_info['title'] = 'Counts distributions'
script_info['script_description'] = 'Dot, line, histogram or box plot '+\
        'of counts distributions of any number of experiments'
script_info['brief_description'] = 'Dot, line, histogram or box plot '+\
        'of counts distributions of any number of experiments'
script_info['usage'] = 'You can exclude genes '+\
        'or force their inclusion by first using gene_overlap.py to '+\
        'generate a gene list and upload it to the DB with '+\
        'add_expression.py'
script_info['version'] = __version__
script_info['authors'] = __author__
script_info['output_description']= 'PDF/PNG/JPEG dot, line, '+\
        'histogram or box plot'
script_info['help_on_no_arguments'] = True

# Process command-line arguments
req_args = ['collections']
opt_args = ['plot_filename', 'title', 'xlabel', 'ylabel', 'fig_height',
        'fig_width', 'plot_type', 'counts_region', 'y_axis_is_log',
        'normalise_by_RPM']
pos_args = []

script_info['args'] = Args(required_args=req_args, optional_args=opt_args,
    positional_args=pos_args)
script_info['required_options'] = script_info['args'].getReqCogentOpts()
script_info['optional_options'] = script_info['args'].getOptCogentOpts()

def make_plot(score_groups, fig_details, plot_type, plot_fn=None):
    """
        Create dot, line, bar plot or histogram of genes in each study by rank
    """
    rr = RunRecord('make_plot')

    fig = pyplot.figure(figsize=(fig_details.x_size, fig_details.y_size))
    ax = fig.add_subplot(111)
    pyplot.title(fig_details.title)
    pyplot.ylabel(fig_details.y_text)
    pyplot.xlabel(fig_details.x_text)

    rcParams['xtick.direction'] = 'out'
    rcParams['ytick.direction'] = 'out'
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    if plot_type == 'dot':
        for scores in score_groups:
            x = [i for i in xrange(len(scores))]
            ax.plot(x, scores, '.')
    elif plot_type == 'line':
        for scores in score_groups:
            x = [i for i in xrange(len(scores))]
            ax.plot(x, scores)
    elif plot_type == 'hist':
        # stacked histogram
        for scores in score_groups:
            pyplot.hist(scores, bins=20)
    elif plot_type == 'box':
        pyplot.boxplot(score_groups)
    else:
        rr.dieOnCritical('plot type not supported', plot_type)

    if plot_fn:
        if plot_fn.lower().endswith('.pdf'):
            pyplot.savefig(plot_fn, format='pdf')
        elif plot_fn.lower().endswith('.jpg') or\
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
        How do counts distributions vary with rank?
    """
    rr = RunRecord('counts_distribution')
    rr.addCommands(sys.argv)
    args = script_info['args'].parse(window_title='Counts Distribution')

    studies = [CountsStudy(fn) for fn in args.collections]

    fig_details = FigureDetails(x_size=args.fig_width, y_size=args.fig_height,
            title=args.title, x_text=args.xlabel, y_text=args.ylabel)

    if args.normalise_by_RPM:
        for study in studies:
            study.normaliseByRPM()

    score_groups = []
    for study in studies:
        score_groups.append(study.scoresAsRankedArray(metric=args.counts_region,
                log2=args.y_axis_is_log))

    make_plot(score_groups, fig_details, args.plot_type, args.plot_filename)

    rr.display()

if __name__ == '__main__':
    main()

