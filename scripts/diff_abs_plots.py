from __future__ import division

import sys
sys.path.extend(['..', '../src'])

from chippy.express.db_query import make_session, get_genes_by_ranked_diff, get_genes_by_ranked_expr
from chippy.util.run_record import RunRecord
from chippy.util.command_args import Args
from matplotlib import pyplot, rcParams

__author__ = 'Cameron Jack'
__copyright__ = 'Copyright 2011-2013, Gavin Huttley, Anuj Pahwa, Cameron Jack'
__credits__ = ['Cameron Jack']
__license__ = 'GPL'
__maintainer__ = 'Cameron Jack'
__email__ = 'cameron.jack@anu.edu.au'
__status__ = 'pre-release'
__version__ = '1.0'

script_info = {}
script_info['script_description'] = 'Creates two plots showing two '+\
        'expression components (x-axis) and their difference (y-axis) '+\
        'as dot-plots'
script_info['brief_description'] = 'dot-plots of absolute expression vs '+\
        'differential expression'
script_info['version'] = __version__
script_info['authors'] = __author__
script_info['output_description']= 'PDF/PNG/JPEG dotplot'

# Process command-line arguments
req_args = ['sample1', 'sample2', 'diff_sample', 'yaxis_units', 'xaxis_units']
opt_args = ['num_genes', 'use_ranks', 'sample_extremes',
        'extremes_colour', 'signif_colour', 'bulk_colour', 'hide_extremes',
        'hide_signif', 'hide_bulk', 'plot1_name', 'plot2_name'
        'xaxis_text1', 'xaxis_text2']
pos_args = ['db_path']

script_info['args'] = Args(required_args=req_args, optional_args=opt_args,
        positional_args=pos_args)
script_info['required_options'] = script_info['args'].req_cogent_opts
script_info['optional_options'] = script_info['args'].opt_cogent_opts

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

    def report_data_counts(self, caller_name):
        """ simple reporting of counts stored """
        rr = RunRecord('report_data_counts')
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

def load_sample_genes(db_path, diff_sample, sample, sample_extremes):
    """ load all portions of diffs into a dict with keys:
        diff_plus1, diff_noSig, diff_minus1, sample_bot,
        sample_mid, sample_top
    """
    rr = RunRecord('load_sample_genes')
    # convert full identifier to stored name
    diff_sample_name = diff_sample.split(' : ')[0]
    sample_name =  sample.split(' : ')[0]

    if sample_extremes > 0.5:
        rr.addWarning('load_sample_genes', 'sample_extremes option '+\
                'must be less than or equal to 0.5', sample_extremes)
        sample_extremes = 1
        rr.addInfo('load_sample_genes', 'setting extremes to default',
                sample_extremes)

    raw_plot_data = RawPlotData(diff_sample_name, sample_name)

    # get diff genes which are significantly up-regulated
    session = make_session(db_path)
    multitest_signif_val = 1
    raw_plot_data.diff_sig_plus1 =\
            get_genes_by_ranked_diff(session, diff_sample_name,
            multitest_signif_val, biotype='protein_coding',
            data_path=None, rank_by='mean', test_run=False)
    session.close()

    # get diff genes which are significantly down-regulated
    session = make_session(db_path)
    multitest_signif_val = -1
    raw_plot_data.diff_sig_minus1 =\
            get_genes_by_ranked_diff(session, diff_sample_name,
            multitest_signif_val, biotype='protein_coding',
            data_path=None, rank_by='mean', test_run=False)
    session.close()

    # get diff genes which are neither up nor down
    session = make_session(db_path)
    multitest_signif_val = 0
    raw_plot_data.diff_sig_zero = get_genes_by_ranked_diff(session,
            diff_sample_name, multitest_signif_val, biotype='protein_coding',
            data_path=None, rank_by='mean', test_run=False)
    session.close()

    # get absolute expression samples
    session = make_session(db_path)
    sample_genes = get_genes_by_ranked_expr(session, sample_name,
        biotype='protein_coding', data_path=None, rank_by='mean',
        test_run=False)
    session.close()
    sample_genes.sort(key=lambda x: x.MeanScore, reverse=True)
    sample_cutoff = int(len(sample_genes) * sample_extremes)
    rr.addInfo('sample cutoff set', sample_cutoff)

    # set absolute expression middle genes
    raw_plot_data.sample_mid =\
            sample_genes[sample_cutoff:len(sample_genes)-sample_cutoff]

    raw_plot_data.sample_top = sample_genes[:sample_cutoff]
    raw_plot_data.sample_bot = sample_genes[-sample_cutoff:]\
            if sample_cutoff else []

    # Report on data sizes loaded
    raw_plot_data.report_data_counts('load_sample_data')

    return raw_plot_data

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

def build_plot_points(raw_plot_data, use_ranks, num_genes):
    """ using common genes only, create x,y,c data points to plot based on
    area restrictions (group_dict), rank or expression values (use_ranks)
    and possibly a limited number of genes (num_genes) """
    rr = RunRecord('build_plot_points')

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
                rr.addWarning('data missing from plot', gene_id)
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
                rr.addWarning('data missing from plot', gene_id)
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
                        colour=c_hex, is_signif=True, is_extreme=True)
            else:
                rr.addWarning('data missing from plot', gene_id)
        else:
            rr.addWarning('data missing from plot', gene_id)

        if plot_dot is not None:
            plot_dot_list.append(plot_dot)
    if num_genes:
        plot_dot_list = plot_dot_list[:num_genes]
    rr.addInfo('number of plot dots', len(plot_dot_list))
    return plot_dot_list

def make_plot(plot_dots, plot_dict, groups_dict):
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



    if plot_dict['out_name'] is not None:
        if '.' in plot_dict['out_name']:
            plot_format = plot_dict['out_name'].split('.')[-1].lower()
        else:
            plot_format = 'png'

        plot_prefix = plot_dict['out_name'].split('.')[:-1]

        plot_fn = plot_prefix + plot_dict['diff_name'] + '_' +\
              plot_dict['sample_name'] + '_diff_abs_expr_plot' +\
              '.' + plot_format

        if plot_prefix == 'pdf':
            pyplot.savefig(plot_fn, image_format='pdf')
        if plot_prefix == 'png':
            pyplot.savefig(plot_fn, image_format='png')
        if plot_prefix == 'jpg' or plot_prefix == '.jpeg':
            pyplot.savefig(plot_fn, image_format='jpg')

    else:
        pyplot.show()


def main():
    """ plot expression scores or ranks in genes of interest for one or 
    more groups"""
    rr = RunRecord('diff_abs_plots')
    rr.addCommands(sys.argv)
    args = script_info['args'].parse(use_scrollbars=True,
            use_save_load_button=True,
            window_title='Difference vs Absolute Expression Plots')

    groups_dict = dict([('extremes_colour', args.extremes_colour),
            ('signif_colour', args.signif_colour),
            ('bulk_colour', args.bulk_colour),
            ('hide_extremes', args.hide_extremes),
            ('hide_signif', args.hide_signif),
            ('hide_bulk', args.hide_bulk)])

    # Should do number restrictions in load step
    # Load all genes into RawPlotData object
    raw_plot_data1 = load_sample_genes(args.db_path, args.diff_sample,
            args.sample1, args.sample_extremes)

    raw_plot_data2 = load_sample_genes(args.db_path, args.diff_sample,
            args.sample2, args.sample_extremes)

    # get back a list of plot_dot objects with 'x', 'y', 'colour', 'area'
    plot_dots1 = build_plot_points(raw_plot_data1, args.ranks,
            args.num_genes)

    plot_dots2 = build_plot_points(raw_plot_data2, args.ranks,
            args.num_genes)

    plot_dict = dict([('out_name', args.plot1_name),
            ('title', args.title),
            ('y_text', args.yaxis_text),
            ('y_units', args.yaxis_units),
            ('x_text', args.xaxis_text),
            ('x_units', args.xaxis_units),
            ('diff_name', raw_plot_data1.diff_name),
            ('sample_name', raw_plot_data1.sample_name)])

    make_plot(plot_dots1, plot_dict, groups_dict)

    plot_dict['sample_name'] = raw_plot_data2.sample_name
    plot_dict['out_name'] = args.plot2_name
    plot_dict['x_axis_text'] = args.xaxis2_text
    plot_dict['x_units'] = args.xaxis2_units

    make_plot(plot_dots2, plot_dict, groups_dict)

    rr.display()

if __name__ == '__main__':
    main()

