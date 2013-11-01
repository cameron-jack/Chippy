from __future__ import division
from math import log10, floor, ceil

import os, sys, glob
sys.path.extend(['..'])

import numpy
import gzip

from optparse import make_option
from chippy.express import db_query
from cogent.util.misc import parse_command_line_parameters
from chippy.util.run_record import RunRecord
from matplotlib import pyplot, rcParams
from chippy.draw.plottable import FigureDetails

__author__ = 'Cameron Jack'
__copyright__ = 'Copyright 2011, Gavin Huttley, Anuj Pahwa, Cameron Jack'
__credits__ = ['Cameron Jack']
__license__ = 'GPL'
__maintainer__ = 'Cameron Jack'
__email__ = 'cameron.jack@anu.edu.au'
__status__ = 'pre-release'
__version__ = '0.2'

def _make_sample_choices(session):
    """returns the available choices for target gene samples"""
    #samples = ['%s : %s' % (s.name, s.description)
    #           for s in db_query.get_target_sample(session)]
    #samples.insert(0, None)
    samples = db_query.get_samples(session)
    if not samples:
        samples = [None]

    return samples

def _create_plot_options_required(session):

    # essential source files
    samples = _make_sample_choices(session)
    opt_sample = make_option('-s', '--sample', type='choice',
                              help='Choose the expression study [default: %default]',
                              choices=[str(s) for s in samples])

    exp_absolute = 'Expression data: absolute ranked'
    exp_diff = 'Expression data: difference in expression between samples'
    target_genes ='Target gene list'
    opt_sample_type = make_option('-y', '--sample_type', type='choice',
        choices=[exp_absolute, exp_diff, target_genes],
        help='Select the type of data you want entered from %s' %\
         str([exp_absolute, exp_diff,target_genes]))
    
    opt_plot_type = make_option('-p', '--plot_type', type='choice',
            choices=['hist', 'box', 'dot'],
            help='Select the type of plot you want entered from %s' %\
             str(['hist', 'box', 'dot']))    

    opt_outfile = make_option('-o', '--outfile', type='string',
        help='Output file path, appends .png automatically')
    
    opt_genefile = make_option('--genefile', type='string', default=None,
        help='Final gene list file output path, as pickle.gz')
    
    opt_orig = make_option('--orig', type='choice',
            help='Choose the original expression study involved in an '\
            'expression_diff experiment [default: %default]',
            choices=[str(s) for s in samples])    

    required_opts = [opt_sample, opt_sample_type, opt_plot_type, opt_outfile,
                     opt_genefile, opt_orig]    

    return required_opts

def _create_extra_sample_options(session):
    # essential source filesx
    samples = _make_sample_choices(session)
    opt_sample = make_option('-s', '--sample', type='choice',
                              help='Choose the expression study [default: %default]',
                              choices=[str(s) for s in samples])
    opt_multitest_signif1 = make_option('-m', '--multitest_signif_val', type='int',
        help='Restrict plot to genes that pass multitest significance,'\
             'valid values: 1, 0, -1', default=None)

    """ ugly hacks for dealing with more than 1 sample """
    # Ugly hack for coping with two samples
    opt_sample2 = make_option('--s2', type='choice',
            help='Choose the expression study [default: %default]',
            choices=[str(s) for s in samples])
    opt_multitest_signif2 = make_option('--m2', type='int',
            help='Restrict plot to genes that pass multitest significance,'\
            'valid values: 1, 0, -1', default=None)
                               
    # Uglier hack for coping with three samples
    opt_sample3 = make_option('--s3', type='choice',
            help='Choose the expression study [default: %default]',
            choices=[str(s) for s in samples])
    opt_multitest_signif3 = make_option('--m3', type='int',
            help='Restrict plot to genes that pass multitest significance,'\
            'valid values: 1, 0, -1', default=None)

    # Even Uglier hack for coping with four samples
    opt_sample3 = make_option('--s4', type='choice',
        help='Choose the expression study [default: %default]',
        choices=[str(s) for s in samples])
    opt_multitest_signif4 = make_option('--m4', type='int',
        help='Restrict plot to genes that pass multitest significance,'\
             'valid values: 1, 0, -1', default=None)
    
    extra_sample_opts = [opt_sample2, opt_sample3, opt_multitest_signif2, opt_multitest_signif3]
    return extra_sample_opts

def _create_sampling_options():
    
    opt_num_genes = make_option('-n', '--num_genes', default=None,
        help='Number of ranked genes to get expression scores for '\
         '[default: %default]')

    opt_ranks = make_option('-r', '--ranks', action='store_true',
        help='Plot expression ranks instead of expression scores', default=False)
    
    sampling_opts = [opt_num_genes, opt_ranks]
    return sampling_opts

def _create_plot_options():
    opt_title = make_option('--title', type='string', default='Info Plot',
        help='Text for the title of the plot [default: %default]')
    opt_yaxis_text = make_option('--yaxis_text', type='string', default='Expression',
        help='Text for y-axis of plot [default: %default]')
    opt_xaxis_text = make_option('--xaxis_text', type='string', default='Expression',
        help='Text for x-axis of plot [default: %default]')

    plot_opts = [opt_title, opt_yaxis_text, opt_xaxis_text]
    return plot_opts

def set_environment():
    """ create the DB session and run options """

    # Create DB session
    if 'CHIPPY_DB' in os.environ:
        db_path = os.environ['CHIPPY_DB']
    else:
        raise RuntimeError('You need to set an environment variable '
                'CHIPPY_DB that indicates where to find the database')
    session = db_query.make_session('sqlite:///%s' % db_path)

    # Describe the application
    script_info = {}
    script_info['title'] = 'Plot gene expression by rank or score'
    script_info['script_description'] = 'Histogram or boxplot of gene '\
            'expression for up to 3 groups of gene expression in the ChipPyDB'
    script_info['usage'] = 'You can exclude genes '\
            'or force their inclusion by first using gene_overlap.py to ' \
            'generate a gene list and upload it to the DB with ' \
            'add_expression.py, then use --iX or --eX to Include or Exclude '\
            'from the corresponding --sX Sample.'
    script_info['version'] = __version__
    script_info['authors'] = __author__
    script_info['output_description']= 'PNG/PDF histogram, boxplot'
    script_info['help_on_no_arguments'] = True

    ### All inputs are divided into logical groupings

    # Required inputs:
    optSet_required_inputs = _create_plot_options_required(session)
    # Extra sample options:
    optSet_extra_inputs = _create_extra_sample_options(session)
    # Sampling is for grouping and filtering by expression:
    optSet_sampling_inputs = _create_sampling_options()    
    # Plot options:
    optSet_plot_inputs = _create_plot_options()

    ### Incorporate all options

    script_info['required_options'] = optSet_required_inputs
    script_info['optional_options'] = optSet_extra_inputs +\
        optSet_sampling_inputs + optSet_plot_inputs

    return db_path, script_info

def load_data(filename):
    if sample_type == 'Expression data: absolute ranked':
        print 'Querying sample'
        sample_name = opts.sample.split(' : ')[0]
        sample_genes = db_query.get_ranked_expression(session, sample_name,
            biotype='protein_coding', data_path=None, rank_by='mean',
            test_run=False)
        rr.addInfo('Expression plots','genes in set 1', len(sample_genes))

        if opts.s2 is not None:
            # IF YOU DON'T DO THIS PYTHON OVERWRITES SAMPLE_GENES! WTF! SERIOUSLY!
            session.close()
            session = db_query.make_session('sqlite:///%s' % db_path)
            print 'Querying sample2'
            sample_name2 = opts.s2.split(' : ')[0]
            sample_genes2 = db_query.get_ranked_expression(session, sample_name2,
                biotype='protein_coding', data_path=None, rank_by='mean',
                test_run=False)
            rr.addInfo('Expression plots','genes in set 2', len(sample_genes2))

    elif sample_type == 'Expression data: difference in expression between samples':
        if (opts.multitest_signif_val is not None) and not\
        (-1 <= opts.multitest_signif_val <= 1):
            raise RuntimeError('multitest_signif_val1 is not -1, 0, 1 or None. '\
                               'Halting execution.')
        if (opts.m2 is not None) and not (-1 <= opts.m2 <= 1):
            raise RuntimeError('multitest_signif_val12is not -1, 0, 1 or None. '\
                               'Halting execution.')

        sample_name = opts.sample.split(' : ')[0]
        print 'Querying sample'
        sample_genes = db_query.get_ranked_expression_diff(session, sample_name,
            opts.multitest_signif_val, biotype='protein_coding',
            data_path=None, rank_by='mean', test_run=False)
        rr.addInfo('Expression plots','diff genes in set 1', len(sample_genes))

        if opts.s2 is not None:
            # IF YOU DON'T DO THIS PYTHON OVERWRITES SAMPLE_GENES! WTF! SERIOUSLY!
            session.close()
            session = db_query.make_session('sqlite:///%s' % db_path)
            sample_name2 = opts.s2.split(' : ')[0]
            print 'Querying sample'
            sample_genes2 = db_query.get_ranked_expression_diff(session, sample_name2,
                opts.m2, biotype='protein_coding',
                data_path=None, rank_by='mean', test_run=False)
            rr.addInfo('Expression plots','diff genes in set 2', len(sample_genes2))
    else:
        print 'Other options not defined yet, choose from %s '\
              'or %s' % (exp_absolute, exp_diff)
        raise RuntimeError ('Incorrect sample type given')

def make_plots(data_list, plot_type, fig_details=None, plot_file=None,
        plot_file_type='png', rr=RunRecord()):

    if not fig_details:
        fig_details = FigureDetails()

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

    elif opts.plot_type == 'box':
        pyplot.boxplot(data_list)

    else:
        rr.addInfo('make_plots', 'plot type not supported', plot_type)
        return rr

    pyplot.show()
    if plot_file:
        outfile = plot_file + plot_type
        pyplot.savefig(outfile, format=plot_type)


def load_all_data(samples, db_path, rr=RunRecord()):
    """ loads both ChipPy counts files and expression records in a ChipPy DB"""

    samples = map(str, samples.strip().split(','))
    sample_names = []
    for sample in samples:
        sample_name = sample.split(' : ')[0]
        sample_names.append(sample_name)

    session = db_query.make_session('sqlite:///%s' % db_path)

    chrm_gene_list = []
    expr_gene_list = []

    for sample in sample_names:
        if os.path.isfile(sample):
            try:
                # to load counts data from file
                file1 = gzip.GzipFile(sample, 'rb')
                data = numpy.load(file1)
                d = data.tolist()
                ranks = d['ranks']
                counts = d['counts']
                labels = d['labels']

                for count, label, rank in zip(counts, labels, ranks):
                    gene_record = ChrmGene(label, sample, counts=count, rank=rank)
                    chrm_gene_list.append(gene_record)
                rr.addInfo('load_all_data', 'genes found in ' + sample,
                        len(labels))

            except IOError: # some exception type
                rr.addError('load_data', 'file found but could not be read', sample)
        else:
            # to load as expression data
            #sample_type == 'Expression data: absolute ranked'
            print 'Querying sample from ChippyDB'
            sample_genes = db_query.get_ranked_expression(session, sample,
                    biotype='protein_coding', data_path=None, rank_by='mean',
                    test_run=False)
            for gene in sample_genes:
                gene_record = ExprGene(gene.stableId, sample, rank=gene.Rank,
                        expr=gene.MeanScore)
                expr_gene_list.append(gene_record)
            rr.addInfo('load_all_data','genes found in ' + sample,
                    len(sample_genes))

    session.close()
    return chrm_gene_list, expr_gene_list, rr

def main():
    """
        Compare the distributions of gene expression relative to
        gene expression rank.
    """
    
    # Get command-line inputs
    db_path, script_info = set_environment()
    option_parser, opts, args =\
                parse_command_line_parameters(**script_info)
    
    rr = RunRecord()

    if opts.sample is None:
        raise RuntimeError('No samples given')

    chrm_gene_list, expr_gene_list, rr = load_all_data(opts.samples, db_path=db_path, rr=rr)

    if opts.num_genes is not None:
        num_genes = int(opts.num_genes)
        
        sample_genes.sort(key = lambda gene: gene.MeanScore)
        if opts.multitest_signif_val != -1:
            sample_genes.reverse()
        sample_genes = sample_genes[:num_genes]
        #if opts.multitest_signif_val == -1:
        #    # need to invert all difference scores
        #    for gene in sample_genes:
        #        gene.Scores = [0 - gene.MeanScore]
        if sample_genes2 is not None:
            sample_genes2.sort(key = lambda gene: gene.MeanScore)
            if opts.m2 != -1:
                sample_genes2.reverse()
            sample_genes2 = sample_genes2[:num_genes]
        

    # Find overlap for set1
    overlapping_gene_ids = set()
    sample_gene_ids1 = []
    sample_gene_ids2 = []
    orig_gene_ids = []
    for gene in sample_genes:
        sample_gene_ids1.append(gene.ensembl_id)
    if sample_genes2 is not None:
        for gene in sample_genes2:
            sample_gene_ids2.append(gene.ensembl_id)
    for gene in orig_genes:
        orig_gene_ids.append(gene.ensembl_id)

    overlapping_gene_ids1 = set(sample_gene_ids1).intersection(set(orig_gene_ids))
    overlapping_gene_ids2 = set(sample_gene_ids2).intersection(set(orig_gene_ids))

    # get matching scores and ids
    sample_scores1 = []
    sample_ids1 = []
    sample_scores2 = []
    sample_ids2 = []
    orig_scores = []
    orig_ids = []
    for gene in sample_genes:
        sample_ids1.append(gene.ensembl_id)
        if opts.ranks:
            sample_scores1.append(gene.Rank)
        else:
            sample_scores1.append(gene.MeanScore)
            
    for gene in sample_genes2:
        sample_ids2.append(gene.ensembl_id)
        if opts.ranks:
            sample_scores2.append(gene.Rank)
        else:
            sample_scores2.append(gene.MeanScore)    

    for gene in orig_genes:
        orig_ids.append(gene.ensembl_id)
        if opts.ranks:
            orig_scores.append(gene.Rank)
        else:
            orig_scores.append(gene.MeanScore)

    sample_dict1 = dict((sample_ids1[j], sample_scores1[j]) for j in range(len(sample_genes)))
    sample_dict2 = dict((sample_ids2[l], sample_scores2[l]) for l in range(len(sample_genes2)))
    orig_dict = dict((orig_ids[k], orig_scores[k]) for k in range(len(orig_genes)))

    x_list = []
    y_list = []
    for id in overlapping_gene_ids1:
        y_list.append(sample_dict1[id])
        x_list.append(orig_dict[id])
        
    for id in overlapping_gene_ids2:
        y_list.append(sample_dict2[id])
        x_list.append(orig_dict[id])

    x_data = numpy.array(x_list)
    y_data = numpy.array(y_list)

    data_list = [x_data, y_data]

    rr= make_plots(data_list, opts.plot_type, rr=rr)

    rr.display()

if __name__ == '__main__':
    main()

