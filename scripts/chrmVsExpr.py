from __future__ import division
from math import log10, floor, ceil

import os, sys, glob, math
sys.path.extend(['..', '../src'])

import numpy
import pickle
import gzip

from optparse import make_option
from chippy.core.collection import RegionCollection
from chippy.express import db_query
from cogent.util.misc import parse_command_line_parameters
from chippy.util.run_record import RunRecord
from matplotlib import pyplot, rcParams
from chippy.core.collection import RegionCollection, column_sum, column_mean, stdev

__author__ = 'Cameron Jack'
__copyright__ = 'Copyright 2011, Gavin Huttley, Anuj Pahwa, Cameron Jack'
__credits__ = ['Gavin Huttley, Cameron Jack']
__license__ = 'GPL'
__maintainer__ = 'Gavin Huttley'
__email__ = 'Gavin.Huttley@anu.edu.au'
__status__ = 'alpha'
__version__ = '0.1'

def _make_sample_choices(session):
    """returns the available choices for target gene samples"""
    #samples = ['%s : %s' % (s.name, s.description)
    #           for s in db_query.get_target_sample(session)]
    #samples.insert(0, None)
    samples = db_query.get_samples(session)
    if not samples:
        samples = [None]

    return [str(s) for s in samples]

def _create_plot_options_required(db_path):
    """ essential sources and options"""
    session = db_query.make_session('sqlite:///%s' % db_path)

    samples = _make_sample_choices(session)
    opt_samples = make_option('-s', '--samples', type='string',
            help='A comma separated list of data files or ChipPyDB ' \
            'experiments. Choose from: ' + '\n'.join([s for s in samples]))
    opt_collections = make_option('-c', '--collections', type='string',
            help='A comma separated list of collection data files.')

    opt_plot_type = make_option('-p', '--plot_type', type='choice',
            choices=['dot', 'line'], help='Select the type of plot you want '\
            'entered from %choices')

    opt_x_axis_type = make_option('-x', '--x_axis_type', type='choice',
            choices=['expression', 'chrm counts', 'expr counts'],
            help='Select the data unit type for the x-axis')

    opt_y_axis_type = make_option('-y', '--y_axis_type', type='choice',
            choices=['expression', 'chrm counts', 'expr counts'],
            help='Select the data unit type for the y-axis')

    opt_output_file = make_option('-o', '--output_file', type='string',
            help='Output plot file path, appends file type automatically')

    required_opts = [opt_samples, opt_collections, opt_plot_type,
            opt_x_axis_type, opt_y_axis_type, opt_output_file]

    return required_opts

def _create_sampling_options(db_path):
    """ options for manipulating interrogated data """
    session = db_query.make_session('sqlite:///%s' % db_path)
    opt_num_genes = make_option('-n', '--num_genes', default=None,
        help='Number of ranked genes to get expression scores for '\
         '[default: %default]')
    opt_bottom_genes = make_option('--bottom', default=False,
            action='store_true', help='Look at bottom n genes rather than '\
            'top genes. [Default: %default]')

    opt_x_ranks = make_option('--x_axis_is_ranks', action='store_true',
            help='Plot x-axis as ranks rather than absolute values',
            default=False)
    opt_y_ranks = make_option('--y_axis_is_ranks', action='store_true',
        help='Plot y-axis as ranks rather than absolute values',
        default=False)

    opt_x_log = make_option('--x_axis_is_log', action='store_true',
            default=False, help='Plot x-axis with log2 values')
    opt_y_log = make_option('--y_axis_is_log', action='store_true',
        default=False, help='Plot y-axis with log2 values')

    samples = _make_sample_choices(session)
    session.close()

    opt_include_genes = make_option('-i', '--include', type='choice',
            help='Choose the target gene list to include [default: %default]',
            choices=[str(s) for s in samples])

    opt_exclude_genes = make_option('-e', '--exclude', type='choice',
            help='Choose the target gene list to exclude [default: %default]',
            choices=samples)
    
    sampling_opts = [opt_num_genes, opt_bottom_genes, opt_x_ranks, opt_y_ranks,
            opt_x_log, opt_y_log, opt_include_genes, opt_exclude_genes]
    return sampling_opts

def _create_plot_options():
    """ options for data display """
    opt_title = make_option('--title', type='string', default='Info Plot',
            help='Text for the title of the plot [default: %default]')
    opt_yaxis_text = make_option('--yaxis_text', type='string', default=None,
            help='Text for y-axis of plot [default: %default]')
    opt_xaxis_text = make_option('--xaxis_text', type='string', default=None,
            help='Text for x-axis of plot [default: %default]')

    opt_fig_height = make_option('-H', '--fig_height', type='float',
            default=2.5*3, help='Figure height (cm) [default: %default]')
    opt_fig_width = make_option('-W', '--fig_width', type='float',
            default=2.5*5, help='Figure width (cm) [default: %default]')

    opt_output_type = make_option('-t', '--output_type', type='choice',
            default = 'png', choices=['png', 'pdf', 'none'],
            help='Plot file format. Choices: %choices [default: %default]')

    plot_opts = [opt_title, opt_yaxis_text, opt_xaxis_text, opt_fig_height,
            opt_fig_width, opt_output_type]
    return plot_opts

def set_environment():
    """ create the DB session and run options """

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

    # link to DB
    if 'CHIPPY_DB' in os.environ:
        db_path = os.environ['CHIPPY_DB']
    else:
        raise RuntimeError('You need to set an environment variable '
                           'CHIPPY_DB that indicates where to find the database')

    ### All inputs are divided into logical groupings

    # Required inputs:
    optSet_required_inputs = _create_plot_options_required(db_path)
    # Sampling is for grouping and filtering by expression:
    optSet_sampling_inputs = _create_sampling_options(db_path)
    # Plot options:
    optSet_plot_inputs = _create_plot_options()

    ### Incorporate all options

    script_info['required_options'] = optSet_required_inputs
    script_info['optional_options'] = optSet_sampling_inputs +\
            optSet_plot_inputs

    return db_path, script_info



class Gene(object):
    """ defined by a stableId in a given study """

    def __init__(self, stableId, study, *args, **kwargs):
        super(Gene, self).__init__(*args, **kwargs)
        self.stableId = stableId
        self.study = study

    def __repr__(self):
        return repr((self.stableId, self.study))

class ChrmGene(Gene):
    """ gene entry from a ChipPy study """
    def __init__(self, counts, *args, **kwargs):
        self.counts = counts # a numpy array
        self.feature_pos = len(self.counts)/2
        self.feature_count = self.counts[self.feature_pos]
        self.total_count = numpy.sum(self.counts)
        self.promoter_counts = numpy.sum(self.counts[:len(counts)/2])
        self.coding_counts = numpy.sum(self.counts[len(counts)/2:])
        super(ChrmGene, self).__init__(*args, **kwargs)

    def __repr__(self):
        return repr((self.counts, self.rank, self.feature_pos,
                self.feature_count))

class ExprGene(Gene):
    """ gene entry from a microarray expression study """
    def __init__(self, expr, *args, **kwargs):
        self.expr = expr
        super(ExprGene, self).__init__(*args, **kwargs)

    def __repr__(self):
        return repr(self.expr, self.rank)

class MatchedStudy(object):
    """ plots of chrm vs expr must have both defined for each plot line """
    def __init__(self, chrm_genes, expr_genes):
        self.chrm_genes = chrm_genes
        self.expr_genes = expr_genes

def load_expr(sample, db_path, rr=None):
    """ loads expression records from a ChippyDB """

    if not rr:
        rr = RunRecord()

    sample_name = sample.split(' : ')[0]

    session = db_query.make_session('sqlite:///%s' % db_path)

    expr_gene_list = []
    #sample_type == 'Expression data: absolute ranked'
    print 'Querying sample from ChippyDB'
    sample_genes = db_query.get_ranked_expression(session, sample_name,
            biotype='protein_coding', data_path=None, rank_by='mean',
            test_run=False)
    for gene in sample_genes:
        gene_record = ExprGene(gene.MeanScore, gene.ensembl_id, sample)
        expr_gene_list.append(gene_record)
    rr.addInfo('load_expression','genes found in ' + sample, len(sample_genes))

    return expr_gene_list, rr

def load_chrm(collection, rr=None):
    """ loads gene entries from a ChipPy collection """

    if not rr:
        rr = RunRecord()

    chrm_gene_list = []
    if os.path.isfile(collection):
        try:
            # to load counts data from file
            file1 = gzip.GzipFile(collection, 'rb')
            data = numpy.load(file1)
            d = data.tolist()
            counts = d['counts']
            labels = d['labels']

            for count, label in zip(counts, labels):
                gene_record = ChrmGene(count, str(label), collection)
                chrm_gene_list.append(gene_record)
            rr.addInfo('load_all_data', 'genes found in ' + collection,
                    len(labels))

        except IOError: # some exception type
            rr.addError('load_data', 'file found but could not be read', collection)
    else:
        rr.addError('load_chrm_counts', 'unrecognised collection file', collection)

    return chrm_gene_list, rr

def load_matched_studies(opts, db_path, rr=None):
    """ each expr sample must be matched to a chrm sample """

    if not rr:
        rr = RunRecord()

    collections = opts.collections.strip().split(',')
    expressions = opts.samples.strip().split(',')
    assert collections > 0, "No collections specified"
    assert expressions > 0, "No expression samples specified"
    assert len(collections) == len(expressions),\
            "number of collections doesn't match number of samples"
    matched_studies = []
    for collection, expression in zip(collections, expressions):
        print collection
        print expression
        chrm_genes, rr = load_chrm(collection, rr)
        expr_genes, rr = load_expr(expression, db_path, rr)
        matched_studies.append(MatchedStudy(chrm_genes, expr_genes))

    return matched_studies, rr

def keep_common_genes(matched_studies, rr=None):
    """ keep only those genes that are common to each study pair """

    if not rr:
        rr = RunRecord()

    # get the intersection of all available stableIds
    final_gene_set = set()
    for matched_study in matched_studies:
        study_gene_set = set()

        for gene in matched_study.expr_genes:
            study_gene_set.add(gene.stableId)
        for gene in matched_study.chrm_genes:
            study_gene_set.add(gene.stableId)
        if not len(final_gene_set):
            for stableId in study_gene_set:
                final_gene_set.add(stableId)
        else:
            final_gene_set = final_gene_set.intersection(study_gene_set)

    # Now keep only those genes that are common to all
    for study in matched_studies:
        for gene in study.expr_genes:
            if not gene.stableId in final_gene_set:
                study.remove(gene)
        for gene in study.chrm_genes:
            if not gene.stableId in final_gene_set:
                study.remove(gene)

    rr.addInfo('keep_common_genes', 'number of kept genes', len(final_gene_set))
    return matched_studies, rr

class Plot_point:
    def __init__(self, x=0, y=0):
        self.x = x
        self.y = y

    def __repr__(self):
        return repr(self.x, self.y)

def create_plot_points(matched_studies, opts, rr=None):
    """ convert values in matched studies to plot points.
    X-axis should always be ordered. We can rely on chrm and expr data
    to not have singletons. Output as list of lists of tuples"""

    if not rr:
        rr = RunRecord()

    plot_point_groups = []
    for matched_study in matched_studies:
        plot_points = {}

        for gene in matched_study.expr_genes:
            if opts.x_axis_type.lower() == 'expression':
                plot_points[gene.stableId] = Plot_point(x=gene.expr)
            else: # opts.x_axis_type.lower() == 'chrm_counts':
                plot_points[gene.stableId] = Plot_point(y=gene.expr)

        for gene in matched_study.chrm_genes:
            if gene.stableId in plot_points:
                if opts.x_axis_type.lower() == 'expression':
                    plot_point = plot_points[gene.stableId]
                    plot_point.y = gene.promoter_counts
                    plot_points[gene.stableId] = plot_point
                else: # opts.x_axis_type.lower() == 'chrm_counts':
                    plot_point = plot_points[gene.stableId]
                    plot_point.x = gene.promoter_counts
                    plot_points[gene.stableId] = plot_point
            else:
                rr.addWarning('create_plot_points',
                        'Singleton in chromatin study', gene.ensembl_id)

        # convert dict to list
        point_list = []
        for key, val in plot_points.iteritems():
            point_list.append((val.x, val.y))

        # sort list based on x values
        point_list = sorted(point_list, key=lambda point: point[0])

        rr.addInfo('create_plot_points', 'number of plot points added', len(point_list))
        plot_points = numpy.array(point_list)
        plot_point_groups.append(plot_points)
    return plot_point_groups, rr

def transform_plot_points(plot_point_groups, opts, rr):
    """ convert axes to ranks or log2 transform """
    # log transform axes if required
    new_point_groups = []
    for point_list in plot_point_groups:
        if opts.x_axis_is_log:
            rr.addInfo('transform_plot_points', 'x-axis transformed', 'log2')
            new_list = []
            for entry in point_list:
                x, y = entry
                if x < 0:
                    rr.addWarning('transform_plot_points', 'negative value pre-log', x)
                x += 1
                new_list.append((math.log(x, 2), y))
            point_list = new_list

        if opts.y_axis_is_log:
            rr.addInfo('transform_plot_points', 'y-axis transformed', 'log2')
            new_list = []
            for entry in point_list:
                x, y = entry
                if y < 0:
                    rr.addWarning('transform_plot_points', 'negative value pre-log', y)
                y += 1
                new_list.append((x, math.log(y, 2)))
            point_list =  new_list

        # rank axes if required
        if opts.x_axis_is_ranks:
            rr.addInfo('transform_plot_points', 'x-axis transformed', 'ranks')
            new_list = []
            for i, entry in enumerate(point_list):
                x, y = entry
                new_list.append((i, y))
            point_list = new_list

        if opts.y_axis_is_ranks:
            rr.addInfo('transform_plot_points', 'y-axis transformed', 'ranks')
            new_list = []
            for i, entry in enumerate(point_list):
                x, y = entry
                new_list.append((x, i))
            point_list =  new_list

        new_point_groups.append(point_list)

    return new_point_groups, rr

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

def make_plots(plot_point_groups, plot_type='dot', plot_file=None,
        output_type=None, fig_details=None, rr=None):
    """ plot_point_groups is a list of 2D numpy arrays """

    if not rr:
        rr = RunRecord()

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

    rr.addInfo('make_plots', 'Output plot', plot_type)

    for study in plot_point_groups:
        x_vals = []
        y_vals = []
        for tup in study:
            x,y = tup
            x_vals.append(x)
            y_vals.append(y)
        if plot_type == 'dot':
            ax.plot(x_vals, y_vals, 'o')
        elif plot_type == 'line':
            ax.plot(x_vals, y_vals)
        else:
            rr.addInfo('make_plots', 'plot type not supported', plot_type)
            return rr

    #if plot_type == 'hist':
    #    # stacked histogram
    #    y_vals_list = data_tuple
    #    pyplot.hist(y_vals_list, bins=10, histtype='barstacked')

    #elif plot_type == 'box':
    #    y_vals_list = data_tuple
    #    pyplot.boxplot(y_vals_list)

    pyplot.show()
    if plot_file and output_type:
        if output_type.lower() != 'none':
            outfile = plot_file + '_' + plot_type + '.' + output_type
            pyplot.savefig(outfile, format=output_type)

    return rr

def main():
    """ Comparative plots of count or rank data for chromatin or expression."""
    
    # Get command-line inputs
    db_path, script_info = set_environment()
    option_parser, opts, args =\
                parse_command_line_parameters(**script_info)
    
    rr = RunRecord()

    # Load all required data
    matched_studies, rr = load_matched_studies(opts, db_path, rr=rr)

    matched_studies, rr = keep_common_genes(matched_studies, rr=rr)

    # create plot points
    plot_point_groups, rr = create_plot_points(matched_studies, opts, rr=rr)

    # transform plot points
    plot_point_groups, rr = transform_plot_points(plot_point_groups, opts, rr=rr)

    # limit to top or bottom n genes
    # Insert a call here to limit displayed genes in some way
    # Right now it's not important so we can add this later
    if opts.num_genes is not None:
        rr.addWarning('main', 'Display limiting', 'Not implemented')

    fig = FigureDetails(x_size=opts.fig_width, y_size=opts.fig_height)

    rr= make_plots(plot_point_groups, plot_type=opts.plot_type,
            plot_file=opts.output_file, output_type=opts.output_type,
            fig_details=fig, rr=rr)

    rr.display()
    table = rr.getMessageTable()
    table.writeToFile('chrmVsExpr_rr.txt', sep='\t')

if __name__ == '__main__':
    main()

