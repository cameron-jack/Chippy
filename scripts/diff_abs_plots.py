from __future__ import division
### Difference vs Absolute dot plots ###
###
### Creates two dot plots of an expression difference set vs its absolute
### expression components.

import os, sys, glob
sys.path.extend(['..', '../src'])

import numpy
import pickle
import gzip


from math import log10, floor, ceil
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
__maintainer__ = 'Cameron Jack'
__email__ = 'cameron.jack@anu.edu.au'
__status__ = 'beta'
__version__ = '1.0'

def _make_sample_choices(session):
    """returns the available choices for external gene samples"""
    samples = ['%s : %s' % (s.name, s.description)
            for s in db_query.get_samples(session)]
    if not samples:
        samples = [None]
    return samples

def _create_plot_options_required(session):
    """ essential source files and other inputs """
    samples = _make_sample_choices(session)
    opt_diff_sample = make_option('-d', '--diff_sample', type='choice',
            help='Choose the expression study [default: %default]',
            choices=[str(s) for s in samples])
    
    opt_abs_sample1 = make_option('-s', '--sample1', type='choice',
            help='Choose the expression study [default: %default]',
            choices=[str(s) for s in samples])

    opt_abs_sample2 = make_option('-t', '--sample2', type='choice',
            help='Choose the expression study [default: %default]',
            choices=[str(s) for s in samples])

    opt_yaxis_units = make_option('--yaxis_units', type='string',
        help='Text showing units of y-axis of plot [default: %default]')
    opt_xaxis_units = make_option('--xaxis_units', type='string',
        help='Text showing units of x-axis of plot [default: %default]')
    opt_xaxis2_units = make_option('--xaxis2_units', type='string',
        help='Text showing units of x-axis of plot2 [default: %default]')

    required_opts = [opt_diff_sample, opt_abs_sample1, opt_abs_sample2,
            opt_yaxis_units, opt_xaxis_units, opt_xaxis2_units]

    return required_opts

def _create_sampling_options():
    
    opt_num_genes = make_option('-n', '--num_genes', type='int', default=None,
            help='Number of ranked genes to get expression scores for '\
            '[default: %default]')

    opt_ranks = make_option('-r', '--use_ranks', action='store_true', default=False,
            help='Plot expression ranks instead of expression scores '\
            '[default: %default]')
    
    opt_sample_extremes = make_option('-e', '--sample_extremes', type='float',
            default=0.0, help='Proportion of least and most absolute ' \
            'expressed genes to treat separately. Set to 0.0 to disable '\
            '[default: %default]')
    
    sampling_opts = [opt_num_genes, opt_ranks, opt_sample_extremes]
    return sampling_opts

def _create_plot_options():
    """ details for the plot """

    # plot text and format
    opt_title = make_option('--title', type='string', default=None,
            help='Text for the title of the plot [default: %default]')
    opt_yaxis_text = make_option('--yaxis_text', type='string', default=None,
            help='Text for y-axis of plot [default: %default]')
    opt_xaxis_text = make_option('--xaxis_text', type='string', default=None,
            help='Text for x-axis of plot [default: %default]')
    opt_xaxis2_text = make_option('--xaxis2_text', type='string', default=None,
        help='Text for x-axis of plot2 [default: %default]')

    opt_format = make_option('--plot_format', type='choice', default='PDF',
            choices=['PNG', 'PDF'],
            help="Select the plot format to output: 'PNG' or "\
            "'PDF' [default: %default]")

    opts_plot_details = [opt_title, opt_yaxis_text, opt_xaxis_text,
            opt_xaxis2_text, opt_format]

    # dot colouring
    opt_extremes_colour = make_option('--extremes_colour', default='blue',
            choices=['blue', 'red', 'yellow', 'green', 'magenta', 'orange',
            'cyan'], help='Colour of dots for absolute expression marked '\
            'as extreme. [default: %default]')
    opt_signif_colour = make_option('--signif_colour', default='blue',
            choices=['blue', 'red', 'yellow', 'green', 'magenta', 'orange',
            'cyan'], help='Colour of dots for difference of expression marked '\
            'as significant. [default: %default]')
    opt_bulk_colour = make_option('--bulk_colour', default='blue',
            choices=['blue', 'red', 'yellow', 'green', 'magenta', 'orange',
            'cyan'], help='Colour of dots for all relatively unexceptional '\
            'expression values. [default: %default]')

    opts_plot_colours = [opt_extremes_colour, opt_signif_colour, opt_bulk_colour]

    # hide unwanted plot areas
    opt_hide_extremes = make_option('--hide_extremes', default=False,
            action='store_true', help='Do not show absolute expression '\
            'considered extreme [default: %default]')
    opt_hide_signif = make_option('--hide_signif', default=False,
            action='store_true', help='Do not show difference expression '\
            'considered significant [default: %default]')
    opt_hide_bulk = make_option('--hide_bulk', default=False,
            action='store_true', help='Do not show expression values'\
            'considered normal [default: %default]')

    opts_plot_hiding = [opt_hide_extremes, opt_hide_signif, opt_hide_bulk]

    plot_opts = opts_plot_details + opts_plot_colours + opts_plot_hiding
    return plot_opts

def _create_extra_options():
    """ misc. options that don't fit other categories """
    opt_genefile = make_option('-g', '--genefile', type='string', default=None,
            help='Annotated gene list file output path, as pickle.gz')
    opt_outfile1 = make_option('-o', '--output_prefix1', type='string',
            default=None, help='Output path prefix for first plot')
    opt_outfile2 = make_option('-p', '--output_prefix2', type='string',
            default=None, help='Output path prefix for second plot')
    
    extra_opts = [opt_genefile, opt_outfile1, opt_outfile2]
    return extra_opts

def _create_session():
    # Create DB session
    if 'CHIPPY_DB' in os.environ:
        db_path = os.environ['CHIPPY_DB']
    else:
        raise RuntimeError('You need to set an environment variable '
                           'CHIPPY_DB that indicates where to find the database')
    session = db_query.make_session('sqlite:///%s' % db_path)

    return session

def set_environment():
    """ create the DB session and run options """

    session = _create_session()

    # Describe the application
    script_info = {}
    script_info['title'] = 'Plot gene expression'
    script_info['script_description'] = 'Creates two dot plots of an '\
            'expression difference set vs its absolute expression components.'
    script_info['version'] = __version__
    script_info['authors'] = __author__
    script_info['output_description']= 'PNG/PDF dotplot'
    script_info['help_on_no_arguments'] = True

    ### All inputs are divided into logical groupings

    # Required inputs:
    optSet_required_inputs = _create_plot_options_required(session)

    # Sampling is for grouping and filtering by expression:
    optSet_sampling = _create_sampling_options()
    # Plot options:
    optSet_plot = _create_plot_options()
    # Extra (misc) options:
    optSet_misc = _create_extra_options()

    ### Incorporate all options

    script_info['required_options'] = optSet_required_inputs
    script_info['optional_options'] = optSet_sampling + optSet_plot +\
                                      optSet_misc

    return session, script_info

class RawPlotData(object):
    """ contains lists of genes from specific plot areas, with diff being on
     the y-axis and absolute on the x-axis"""

    diff_name = None
    diff_sig_plus1 = [] # top of plot
    diff_sig_zero = [] # middle of plot
    diff_sig_minus1 = [] # bottom of plot

    sample_name = None
    sample_top = [] # right of plot
    sample_mid = [] # middle of plot
    sample_bot = [] # left of plot

    def __init__(self, diff_name, sample_name):
        super(RawPlotData, self).__init__()
        self.diff_name = diff_name
        self.sample_name = sample_name

    def get_diff_sig_plus1_geneID_set(self):
        diff_sig_plus1_ids = set()
        for gene in self.diff_sig_plus1:
            diff_sig_plus1_ids.add(gene.ensembl_id)
        return diff_sig_plus1_ids

    def get_diff_sig_zero_geneID_set(self):
        diff_sig_zero_ids = set()
        for gene in self.diff_sig_zero:
            diff_sig_zero_ids.add(gene.ensembl_id)
        return diff_sig_zero_ids

    def get_diff_sig_minus1_geneID_set(self):
        diff_sig_minus1_ids = set()
        for gene in self.diff_sig_minus1:
            diff_sig_minus1_ids.add(gene.ensembl_id)
        return diff_sig_minus1_ids

    def get_sample_top_geneID_set(self):
        sample_top_ids = set()
        for gene in self.sample_top:
            sample_top_ids.add(gene.ensembl_id)
        return sample_top_ids

    def get_sample_mid_geneID_set(self):
        sample_mid_ids = set()
        for gene in self.sample_mid:
            sample_mid_ids.add(gene.ensembl_id)
        return sample_mid_ids

    def get_sample_bot_geneID_set(self):
        sample_bot_ids = set()
        for gene in self.sample_bot:
            sample_bot_ids.add(gene.ensembl_id)
        return sample_bot_ids

    def get_diff_geneID_set(self):
        diff_gene_ids = self.get_diff_sig_plus1_geneID_set()
        diff_gene_ids = diff_gene_ids.union(self.get_diff_sig_zero_geneID_set())
        diff_gene_ids = diff_gene_ids.union(self.get_diff_sig_minus1_geneID_set())
        return diff_gene_ids

    def get_sample_geneID_set(self):
        sample_gene_ids = self.get_sample_top_geneID_set()
        sample_gene_ids = sample_gene_ids.union(self.get_sample_mid_geneID_set())
        sample_gene_ids = sample_gene_ids.union(self.get_sample_bot_geneID_set())
        return sample_gene_ids


    def report_data_counts(self, caller_name, rr):
        """ simple reporting of counts stored """
        # Report diff counts
        rr.addInfo(caller_name, 'Difference sample name', self.diff_name)
        rr.addInfo(caller_name, 'diff genes for signif 1',
                len(self.diff_sig_plus1))
        rr.addInfo(caller_name, 'diff genes for signif 0',
                len(self.diff_sig_zero))
        rr.addInfo(caller_name, 'diff genes for signif -1',
                len(self.diff_sig_minus1))

        # Report sample1 counts
        rr.addInfo(caller_name, 'Absolute sample name', self.sample_name)
        rr.addInfo(caller_name, 'top extreme genes for sample 1',
                len(self.sample_top))
        rr.addInfo(caller_name, 'bulk, non-extreme genes for sample 1',
                len(self.sample_mid))
        rr.addInfo(caller_name, 'bottom extreme genes for sample 1',
                len(self.sample_bot))

        return rr

def load_sample_genes(session, diff_sample, sample, sample_extremes,
        groups_dict, rr):
    """ load all portions of diffs into a dict with keys:
    diff_plus1, diff_noSig, diff_minus1, sample_bot, sample_mid, sample_top """

    # convert full identifier to stored name
    diff_sample_name = diff_sample.split(' : ')[0]
    sample_name =  sample.split(' : ')[0]

    if sample_extremes > 0.5:
        rr.addWarning('load_sample_genes', 'sample_extremes option '\
                'must be less than or equal to 0.5', sample_extremes)
        sample_extremes = 1
        rr.addInfo('load_sample_genes', 'setting extremes to default',
                sample_extremes)

    raw_plot_data = RawPlotData(diff_sample_name, sample_name)

    # get diff genes with signif +1
    multitest_signif_val = 1
    raw_plot_data.diff_sig_plus1 = db_query.get_ranked_expression_diff(session,
            diff_sample_name, multitest_signif_val, biotype='protein_coding',
            data_path=None, rank_by='mean', test_run=False)

    # get diff genes with signif -1
    multitest_signif_val = -1
    raw_plot_data.diff_sig_minus1 = db_query.get_ranked_expression_diff(session,
            diff_sample_name, multitest_signif_val, biotype='protein_coding',
            data_path=None, rank_by='mean', test_run=False)

    # get diff genes with signif 0
    multitest_signif_val = 0
    raw_plot_data.diff_sig_zero = db_query.get_ranked_expression_diff(session,
            diff_sample_name, multitest_signif_val, biotype='protein_coding',
            data_path=None, rank_by='mean', test_run=False)

    session.close()
    # get absolute expression samples
    session = _create_session()
    sample_genes = db_query.get_ranked_expression(session, sample_name,
        biotype='protein_coding', data_path=None, rank_by='mean',
        test_run=False)
    sample_genes.sort(key=lambda x: x.MeanScore, reverse=True)
    sample_cutoff = int(len(sample_genes) * sample_extremes)
    #print 'sample cutoff is: %d' % sample_cutoff

    # set absolute expression middle genes
    raw_plot_data.sample_mid = sample_genes[sample_cutoff:len(sample_genes)-sample_cutoff]

    raw_plot_data.sample_top = sample_genes[:sample_cutoff]
    raw_plot_data.sample_bot = sample_genes[-sample_cutoff:] if sample_cutoff else []

    # Report on data sizes loaded
    raw_plot_data.report_data_counts('load_sample_data', rr)

    return raw_plot_data, rr

class PlotDot(object):
    """ defines a data point to be used on a dot 2D plot """

    def __init__(self, x=None, y=None, colour=None, is_signif=False, is_extreme=False):
        super(PlotDot, self).__init__()
        self.x = x
        self.y = y
        self.colour = colour
        self.is_signif = is_signif
        self.is_extreme = is_extreme

def _get_ids_scores(data_area, use_ranks):
    """ returns a dict of id and score """
    ids_score_dict = {}

    for gene in data_area:
        if use_ranks:
            ids_score_dict[gene.ensembl_id] = gene.Rank
        else:
            ids_score_dict[gene.ensembl_id] = gene.MeanScore

    return ids_score_dict

def build_plot_points(raw_plot_data, use_ranks, num_genes, rr):
    """ using common genes only, create x,y,c data points to plot based on
    area restrictions (group_dict), rank or expression values (use_ranks)
    and possibly a limited number of genes (num_genes) """

    # Firstly get ids and matching scores for diff and sample

    diff_plus1_ids_scores = \
            _get_ids_scores(raw_plot_data.diff_sig_plus1, use_ranks)
    diff_zero_ids_scores =\
            _get_ids_scores(raw_plot_data.diff_sig_zero, use_ranks)
    diff_minus1_ids_scores =\
            _get_ids_scores(raw_plot_data.diff_sig_minus1, use_ranks)

    sample_top_ids_scores = \
            _get_ids_scores(raw_plot_data.sample_top, use_ranks)
    sample_mid_ids_scores =\
            _get_ids_scores(raw_plot_data.sample_mid, use_ranks)
    sample_bot_ids_scores =\
            _get_ids_scores(raw_plot_data.sample_bot, use_ranks)

    # Now find intersection of genes for X and Y axes
    diff_gene_set = raw_plot_data.get_diff_geneID_set()
    sample_gene_set = raw_plot_data.get_sample_geneID_set()

    #print 'number of genes in diff set %d' % len(diff_gene_set)
    #print 'number of genes in sample set %d' % len(sample_gene_set)

    diff_sample_common_genes = diff_gene_set.intersection(sample_gene_set)

    # Then use the intersection to pull X/Y scores from the diff and sample
    # dicts by id.

    plot_dot_list =[]

    #print 'number of common genes is %d' % len(diff_sample_common_genes)

    diff_sig_plus1_geneID_set = raw_plot_data.get_diff_sig_plus1_geneID_set()
    diff_sig_zero_geneID_set = raw_plot_data.get_diff_sig_zero_geneID_set()
    diff_sig_minus1_geneID_set = raw_plot_data.get_diff_sig_minus1_geneID_set()
    sample_top_geneID_set = raw_plot_data.get_sample_top_geneID_set()
    sample_mid_geneID_set = raw_plot_data.get_sample_mid_geneID_set()
    sample_bot_geneID_set = raw_plot_data.get_sample_bot_geneID_set()

    for gene_id in diff_sample_common_genes:
        plot_dot = None
        if gene_id in diff_sig_plus1_geneID_set:
            # top right hand of plot
            if gene_id in sample_top_geneID_set:
                c_hex = '#%02x%02x%02x' % (255, 0, 255)
                plot_dot = PlotDot(x=sample_top_ids_scores[gene_id],
                        y=diff_plus1_ids_scores[gene_id],
                        colour=c_hex, is_signif=True, is_extreme=True)
            # top middle of plot
            elif gene_id in sample_mid_geneID_set:
                c_hex = '#%02x%02x%02x' % (0, 0, 255)
                plot_dot = PlotDot(x=sample_mid_ids_scores[gene_id],
                        y=diff_plus1_ids_scores[gene_id],
                        colour=c_hex, is_signif=True)
            # top left hand of plot
            elif gene_id in sample_bot_geneID_set:
                c_hex = '#%02x%02x%02x' % (0, 255, 255)
                plot_dot = PlotDot(x=sample_bot_ids_scores[gene_id],
                        y=diff_plus1_ids_scores[gene_id],
                        colour=c_hex, is_signif=True, is_extreme=True)
            else:
                rr.addWarning('build_plot_points', 'data missing from plot',
                        gene_id)
        elif gene_id in diff_sig_zero_geneID_set:
            # middle, right of plot
            if gene_id in sample_top_geneID_set:
                c_hex = '#%02x%02x%02x' % (255, 0, 0)
                plot_dot = PlotDot(x=sample_top_ids_scores[gene_id],
                        y=diff_zero_ids_scores[gene_id],
                        colour=c_hex, is_extreme=True)
            # middle of plot
            elif gene_id in sample_mid_geneID_set:
                c_hex = '#%02x%02x%02x' % (127,127,127)
                plot_dot = PlotDot(x=sample_mid_ids_scores[gene_id],
                        y=diff_zero_ids_scores[gene_id],
                        colour=c_hex)
            # middle, left of plot
            elif gene_id in sample_bot_geneID_set:
                c_hex = '#%02x%02x%02x' % (0, 255, 0)
                plot_dot = PlotDot(x=sample_bot_ids_scores[gene_id],
                        y=diff_zero_ids_scores[gene_id],
                        colour=c_hex, is_extreme=True)
            else:
                rr.addWarning('build_plot_points', 'data missing from plot',
                    gene_id)
        elif gene_id in diff_sig_minus1_geneID_set:
            # bottom, right of plot
            if gene_id in sample_top_geneID_set:
                c_hex = '#%02x%02x%02x' % (255, 0, 255)
                plot_dot = PlotDot(x=sample_top_ids_scores[gene_id],
                        y=diff_minus1_ids_scores[gene_id],
                        colour=c_hex, is_signif=True, is_extreme=True)
            # bottom, middle
            elif gene_id in sample_mid_geneID_set:
                c_hex = '#%02x%02x%02x' % (0, 0, 255)
                plot_dot = PlotDot(x=sample_mid_ids_scores[gene_id],
                        y=diff_minus1_ids_scores[gene_id],
                        colour=c_hex, is_signif=True)
            # bottom, left
            elif gene_id in sample_bot_geneID_set:
                c_hex = '#%02x%02x%02x' % (0, 255, 255)
                plot_dot = PlotDot(x=sample_bot_ids_scores[gene_id],
                        y=diff_minus1_ids_scores[gene_id],
                        colour='c_hex', is_signif=True, is_extreme=True)
            else:
                rr.addWarning('build_plot_points', 'data missing from plot',
                    gene_id)
        else:
            rr.addWarning('build_plot_points', 'data missing from plot',
                gene_id)

        if plot_dot is not None:
            plot_dot_list.append(plot_dot)

    return plot_dot_list, rr

def make_plot(plot_dots, plot_dict, groups_dict, rr):
    """ generate the plot and save to file """

    # Plotting with same look as Rohan's diff-expression plots
    rcParams['xtick.direction'] = 'out'
    rcParams['ytick.direction'] = 'out'
    fig = pyplot.figure(figsize=(6,8))
    ax = fig.add_subplot(111)

    hide_bulk = groups_dict['hide_bulk']
    hide_signif = groups_dict['hide_signif']
    hide_extremes = groups_dict['hide_extremes']

    for dot in plot_dots:
        if (not hide_bulk \
                or (hide_bulk and (dot.is_signif or dot.is_extreme))) \
                and (not hide_signif \
                or (hide_signif and not dot.is_signif)) \
                and ((not hide_extremes) \
                or (hide_extremes and not dot.is_extreme)):
            ax.scatter(dot.x, dot.y, c=dot.colour, marker='o')


    #pyplot.title(plot_dict['title'])
    if plot_dict['y_text'] is None:
        y_axis_text = plot_dict['diff_name'] + ' ' + plot_dict['y_units']
    else:
        y_axis_text = plot_dict['y_text'] + ' ' + plot_dict['y_units']

    if plot_dict['x_text'] is None:
        x_axis_text = plot_dict['sample_name'] + ' ' + plot_dict['x_units']
    else:
        x_axis_text = plot_dict['x_text'] + ' ' + plot_dict['x_units']

    pyplot.ylabel(y_axis_text)
    pyplot.xlabel(x_axis_text)

    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    pyplot.show()

    if plot_dict['prefix'] is None:
        prefix = ''
    else:
        prefix = [plot_dict] + '_'

    outfile = prefix + plot_dict['diff_name'] + '_' + \
              plot_dict['sample_name'] + '_diff_abs_expr_plot' + \
              '.' + plot_dict['format'].lower()

    pyplot.savefig(outfile, format=plot_dict['format'])

    return rr

def main():
    """ plot expression scores or ranks in genes of interest for one or 
    more groups"""
    
    # Get command-line inputs
    session, script_info = set_environment()
    option_parser, opts, args =\
                parse_command_line_parameters(**script_info)
    
    rr = RunRecord()

    groups_dict = dict([('extremes_colour', opts.extremes_colour),
            ('signif_colour', opts.signif_colour),
            ('bulk_colour', opts.bulk_colour),
            ('hide_extremes', opts.hide_extremes),
            ('hide_signif', opts.hide_signif),
            ('hide_bulk', opts.hide_bulk)])

    # Should do number restrictions in load step
    # Load all genes into RawPlotData object
    raw_plot_data1, rr = load_sample_genes(session, opts.diff_sample,
            opts.sample1, opts.sample_extremes, groups_dict, rr)

    raw_plot_data2, rr = load_sample_genes(session, opts.diff_sample,
            opts.sample2, opts.sample_extremes, groups_dict, rr)

    session.close() # End DB session now data had been loaded

    # get back a list of plot_dot objects with 'x', 'y', 'colour', 'area'
    plot_dots1, rr = build_plot_points(raw_plot_data1, opts.use_ranks,
            opts.num_genes, rr)

    plot_dots2, rr = build_plot_points(raw_plot_data2, opts.use_ranks,
            opts.num_genes, rr)

    plot_dict = dict([('prefix', opts.output_prefix1),
            ('format', opts.plot_format),
            ('title', opts.title),
            ('y_text', opts.yaxis_text),
            ('y_units', opts.yaxis_units),
            ('x_text', opts.xaxis_text),
            ('x_units', opts.xaxis_units),
            ('diff_name', raw_plot_data1.diff_name),
            ('sample_name', raw_plot_data1.sample_name)])

    make_plot(plot_dots1, plot_dict, groups_dict, rr)

    plot_dict['sample_name'] = raw_plot_data2.sample_name
    plot_dict['prefix'] = opts.output_prefix2
    plot_dict['x_axis_text'] = opts.xaxis2_text
    plot_dict['x_units'] = opts.xaxis2_units

    make_plot(plot_dots2, plot_dict, groups_dict, rr)

    rr.display()

if __name__ == '__main__':
    main()

