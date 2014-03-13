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
opt_args = ['plot_name', 'xlabel', 'title', 'ylabel']
pos_args = ['db_path']

script_info['args'] = Args(required_args=req_args, optional_args=opt_args,
    positional_args=pos_args)
script_info['required_options'] = script_info['args'].getReqCogentargs()
script_info['optional_options'] = script_info['args'].getOptCogentargs()

def make_plot(studies, args):
    """
        Create dot, line, bar plot or histogram of genes in each study by rank
    """
    rr = RunRecord('make_plot')

    fig = pyplot.figure(figsize=(fig_details.x_size, fig_details.y_size))
    ax = fig.add_subplot(111)
    pyplot.title(fig_details.title)
    pyplot.ylabel(fig_details.yaxis_text)
    pyplot.xlabel(fig_details.xaxis_text)

    rr.addInfo('top_gene_info', 'Output plot', plot_type)

    if plot_type == 'dot' and len(data_list) == 2:
        rcParams['xtick.direction'] = 'out'
        rcParams['ytick.direction'] = 'out'
        ax.plot(data_list[0], data_list[1], 'o')
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')

    if plot_type == 'hist':
        # stacked histogram
        pyplot.hist(data_list, bins=10, histtype='barstacked')

    elif args.plot_type == 'box':
        pyplot.boxplot(data_list)

    else:
        rr.addInfo('make_plots', 'plot type not supported', plot_type)

    if plot_file:
        outfile = plot_file + plot_type
        pyplot.savefig(outfile, format=plot_type)
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

    make_plot(studies, args)

    rr.display()

if __name__ == '__main__':
    main()

