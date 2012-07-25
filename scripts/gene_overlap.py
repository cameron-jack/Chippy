from __future__ import division
from math import log10, floor, ceil

### Gene_overlap.py
#
# This was used to find the overlapping genes from diff studies that were neither
# particularly up or down regulated across cell cycles. The top 100 of these were
# taken to be examples of actively expressed house-keeping genes.

import os, sys, glob
sys.path.extend(['..', '../src'])

import numpy
import pickle
import gzip

from optparse import make_option
from chippy.core.collection import RegionCollection
from chippy.express import db_query
from cogent.util.misc import parse_command_line_parameters
from chippy.util.run_record import RunRecord
from matplotlib import pyplot

__author__ = 'Cameron Jack'
__copyright__ = 'Copyright 2011, Gavin Huttley, Cameron Jack, Anuj Pahwa'
__credits__ = ['Gavin Huttley, Cameron Jack']
__license__ = 'GPL'
__maintainer__ = 'Cameron Jack'
__email__ = 'Cameron.Jack@anu.edu.au'
__status__ = 'alpha'
__version__ = '0.1'

def get_sample_name(sample):
    """returns sample name from a 'sample : description' string"""
    if str(sample) != 'None':
        sample = sample.split(':')[0].strip()
    else:
        sample = None
    return sample

def _make_sample_choices(session):
    """returns the available choices for target gene samples"""
    samples = ['%s : %s' % (s.name, s.description)
               for s in db_query.get_samples(session)]
    if not samples:
        samples = [None]
    return samples

def _create_required(session):
    """ essential source files and other inputs """
    samples = _make_sample_choices(session)

    opt_sample1 = make_option('-s', '--sample1', type='choice',
            help='Choose the expression study [default: %default]',
            choices=[str(s) for s in samples])
    opt_sample2 = make_option('-t', '--sample2', type='choice',
            help='Choose the expression study [default: %default]',
            choices=[str(s) for s in samples])

    exp_absolute = 'Expression data: absolute ranked'
    exp_diff = 'Expression data: difference in expression between samples'
    target_genes ='Target gene list'
    opt_sample1_type = make_option('-w', '--sample1_type', type='choice',
            choices=[exp_absolute, exp_diff, target_genes],
            help='Select the type of data you want entered from %s' %\
            str([exp_absolute, exp_diff, target_genes]))

    opt_sample2_type = make_option('-x', '--sample2_type', type='choice',
            choices=[exp_absolute, exp_diff, target_genes],
            help='Select the type of data you want entered from %s' %\
            str([exp_absolute, exp_diff, target_genes]))

    intersection_text = 'Intersection: the genes in common between samples'
    union_text = 'Union: the superset of all genes found in given samples'
    complement_text = 'Complement: all genes NOT in common between all samples'
    specific_text = 'Specific: genes that are expressed in only one sample'
    opt_comparison_type = make_option('-c', '--comparison_type', type='choice',
            choices=[intersection_text, union_text, complement_text,
                     specific_text],
            help='Select the type of comparison you want to conduct from %s' %\
            str([intersection_text, union_text, complement_text,
                    specific_text]))

    opt_genefile = make_option('--genefile', type='string',
            help='Final gene list file output path. Text file with one '\
            'stableID per line')

    required_opts = [opt_sample1, opt_sample2, opt_sample1_type,
            opt_sample2_type, opt_comparison_type, opt_genefile]

    return required_opts

def _create_extra_input_options(session):
    samples = _make_sample_choices(session)
    opt_sample3 = make_option('-u', '--sample3', type='choice',
            help='Choose the expression study [default: %default]',
            choices=[str(s) for s in samples])

    exp_absolute = 'Expression data: absolute ranked'
    exp_diff = 'Expression data: difference in expression between samples'
    target_genes ='Target gene list'
    opt_sample3_type = make_option('-y', '--sample3_type', type='choice',
            choices=[exp_absolute, exp_diff, target_genes],
            help='Select the type of data you want entered from %s' %\
            str([exp_absolute, exp_diff, target_genes]))

    opt_expression_sample1 = make_option('--expression_sample1',
            type='choice', default=None, help='Choose the expression study '\
            'matching sample1 (for when you want to select by top expressing ' \
            'genes for instance [default: %default]',
            choices=[str(s) for s in samples])
    opt_expression_sample2 = make_option('--expression_sample2',
            type='choice', default=None, help='Choose the expression study ' \
            'matching sample2 (for when you want to select by top expressing ' \
            'genes for instance [default: %default]',
            choices=[str(s) for s in samples])
    opt_expression_sample3 = make_option('--expression_sample3',
            type='choice', default=None, help='Choose the expression study ' \
            'matching sample3 (for when you want to select by top expressing ' \
            'genes for instance [default: %default]',
            choices=[str(s) for s in samples])

    opt_favoured_expression_sample = make_option('--favoured_expression_sample',
            type='int', default=1, help='Whenever a gene in found in multiple' \
            'studies, choose which numbered expression study to draw expressed' \
            'rank from. [default: %default]')

    extra_input_opts = [opt_sample3, opt_sample3_type, opt_expression_sample1,
            opt_expression_sample2, opt_expression_sample3,
            opt_favoured_expression_sample]

    return extra_input_opts

def _create_sampling_options():
    """ these are options that restrict the sample sets in various ways """

    opt_multitest_signif1 = make_option('--m1', type='int',
        help='Restrict plot to genes that pass multitest significance,'\
             'valid values: 1, 0, -1', default=None)
    opt_multitest_signif2 = make_option('--m2', type='int',
        help='Restrict plot to genes that pass multitest significance,'\
             'valid values: 1, 0, -1', default=None)
    opt_multitest_signif3 = make_option('--m3', type='int',
        help='Restrict plot to genes that pass multitest significance,'\
             'valid values: 1, 0, -1', default=None)

    opt_num_genes = make_option('-n', '--num_genes', type='int', default=None,
            help='Number of ranked genes to get expression scores for. You ' \
            'must also give --expression_sampleX for each sample so that ' \
            'expression scores can be selected. [default: %default]')

    opt_sample_extremes = make_option('-e', '--sample_extremes', type='float',
            default=0.0, help='Proportion of least and most absolute '\
            'expressed genes to treat separately. Set to 0.0 to disable '\
            '[default: %default]')
    opts_top_extreme_off = make_option('--top_extreme_off',
            action='store_true', default='False', help='If you set sample '\
            'extremes then by default both extremes are kept and bulk '\
            'expressing genes are dropped. This will disable the high '\
            'expressing portion of extreme expressing genes')
    opts_bottom_extreme_off = make_option('--bottom_extreme_off',
            action='store_true', default='False', help='If you set sample '\
            'extremes then by default both extremes are kept and bulk '\
            'expressing genes are dropped. This will disable the low '\
            'expressing portion of extreme expressing genes')

    sampling_opts = [opt_num_genes, opt_sample_extremes, opt_multitest_signif1,
                     opt_multitest_signif2, opt_multitest_signif3,
                     opts_top_extreme_off, opts_bottom_extreme_off]
    return sampling_opts

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
    script_info['title'] = 'Gene overlap exporter'
    script_info['script_description'] = 'Investigate intersections or unions '\
            'between up to 3 expression or expression_diff databases by rank '\
            'or true expression.'
    script_info['version'] = __version__
    script_info['authors'] = __author__
    script_info['output_description']= 'Zip and pickled list of genes'
    script_info['help_on_no_arguments'] = True

    ### All inputs are divided into logical groupings

    # Required inputs:
    optSet_required_inputs = _create_required(session)

    # Options for additional inputs
    optSet_extra_inputs = _create_extra_input_options(session)

    # Sampling is for grouping and filtering by expression:
    optSet_sampling_inputs = _create_sampling_options()

    ### Incorporate all options

    script_info['required_options'] = optSet_required_inputs
    script_info['optional_options'] = optSet_extra_inputs + optSet_sampling_inputs

    return script_info

def getExpressedGenes(session, sample, sample_type='Expression data: '\
        'absolute ranked', multitest_signif_val=None, sample_extremes=0.0,
        top_extreme_off=False, bottom_extreme_off=False, rr=RunRecord()):
    """ Return stableIDs for genes of interest """
    sample_name = sample.split(' : ')[0]
    #print sample_name + '\n'

    if sample_type == 'Expression data: absolute ranked':
        sample_genes = db_query.get_ranked_expression(session, sample_name,
                biotype='protein_coding', data_path=None, rank_by='mean',
                test_run=False)

        if sample_extremes > 0.0:
            sample_genes.sort(key=lambda x: x.MeanScore, reverse=True)
            sample_cutoff = int(len(sample_genes) * sample_extremes)
            rr.addInfo('getExpressionIDs' ,'sample cutoff', sample_cutoff)

            # set absolute expression gene regions
            # sample_mid = sample_genes[sample_cutoff:len(sample_genes)-sample_cutoff]
            if not top_extreme_off:
                sample_top = sample_genes[:sample_cutoff]
            else:
                sample_top = []

            if not bottom_extreme_off:
                sample_bottom = sample_genes[-sample_cutoff:] \
                        if sample_cutoff else []
            else:
                sample_bottom = []

            sample_genes = sample_top + sample_bottom

    elif sample_type == 'Expression data: difference in expression between samples':
        if (multitest_signif_val is not None) and not \
                (-1 <= multitest_signif_val <= 1):
            raise RuntimeError('multitest_signif_val1 is not -1, 0, 1 or None. '\
                    'Halting execution.')

        sample_genes = db_query.get_ranked_expression_diff(session, sample_name,
                multitest_signif_val, biotype='protein_coding',
                data_path=None, rank_by='mean', test_run=False)

    else:
        raise RuntimeError('Target genes not implemented yet, or %d') % sample_type

    if not len(sample_genes):
        rr.addError('getExpressionIDs', 'Sample genes remaining', 0)
        rr.display()
        raise RuntimeError('No sample genes in final set')

    rr.addInfo('getExpressedGenes', sample_name + ' genes returned', len(sample_genes))

    return sample_genes, rr

def compare_ids(sample1_ids, sample2_ids, sample3_ids, comparison_type, rr):
    if comparison_type == 'Intersection':
        output_id_set = sample1_ids.intersection(sample2_ids)
        if len(sample3_ids) > 0:
            output_id_set = output_id_set.intersection(sample3_ids)

    elif comparison_type == 'Union':
        output_id_set = sample1_ids.union(sample2_ids)
        output_id_set = output_id_set.union(sample3_ids)

    elif comparison_type == 'Complement':
        # Those genes that do not intersect in all samples
        unionSet = sample1_ids.union(sample2_ids.union(sample3_ids))
        intersectionSet = sample1_ids.intersection(sample2_ids)
        if len(sample3_ids) > 0:
            intersectionSet = intersectionSet.intersection(sample3_ids)
        output_id_set = unionSet - intersectionSet

    elif comparison_type == 'Specific':
        # Those genes found in only one study
        if not len(sample3_ids):
            output_id_set = sample1_ids.union(sample2_ids) -\
                            sample1_ids.intersection(sample2_ids)
        else:
            output_id_set = sample1_ids - sample2_ids - sample3_ids
            output_id_set = output_id_set.union(sample2_ids - sample1_ids -\
                                                sample3_ids)
            output_id_set = output_id_set.union(sample3_ids - sample1_ids -\
                                                sample2_ids)

    else:
        rr.display()
        raise RuntimeError('Invalid comparison type %d') % comparison_type

    return output_id_set, rr

def truncate_ids_by_expressed(output_id_set, expressed_gene_list,
        num_genes, rr):
    """ Return output_id_set ranked and truncated by expression data """

    # create a dict of id and score so that duplicates are removed
    id_score_dict = {}
    for gene in expressed_gene_list:
        if gene.ensembl_id in output_id_set:
            id_score_dict[gene.ensembl_id] = gene.MeanScore

    # Now we have a single set of ids, we need to get a list of IDs ordered
            # by expression
    ordered_id_list = []
    for geneID, score in sorted(id_score_dict.iteritems(),
            key=lambda (id, score): (score, id), reverse=True):
        ordered_id_list.append(geneID)

    truncated_geneID_list = ordered_id_list[0:num_genes]

    rr.addInfo('truncate_ids_by_expressed', 'geneID list truncated to',
            len(truncated_geneID_list))

    return_id_set = set(truncated_geneID_list)

    return return_id_set, rr

def restrict_by_num_genes(output_id_set, num_genes, expression_sample1,
        expression_sample2, expression_sample3, favoured_expression, rr):
    """ rank output_id_set by expression values and truncate """

    if not (expression_sample1 and expression_sample2) or \
            (favoured_expression == 3 and not expression_sample3):
        rr.display()
        raise RuntimeError('Expression samples required to match each input')

    # Load all expressed gene lists
    session = _create_session()
    expression_genes1, rr = getExpressedGenes(session, expression_sample1, rr=rr)
    session.close()
    session = _create_session()
    expression_genes2, rr = getExpressedGenes(session, expression_sample2, rr=rr)
    session.close()
    expression_genes3 = []
    if expression_sample3:
        session = _create_session()
        expression_genes3, rr = getExpressedGenes(session, expression_sample3, rr=rr)
        session.close()

    # order total gene list by favoured expression
    if favoured_expression == 1:
        expressed_gene_list = expression_genes2 + expression_genes3 + \
                              expression_genes1
    elif favoured_expression == 2:
        expressed_gene_list = expression_genes1 + expression_genes3 + \
                              expression_genes2
    elif favoured_expression == 3:
        expressed_gene_list = expression_genes1 + expression_genes2 + \
                              expression_genes3

    rr.addInfo('restrict_by_num_genes', 'total expression values loaded',
            len(expressed_gene_list))

    # Rank and truncate
    output_id_set, rr = truncate_ids_by_expressed(output_id_set,
            expressed_gene_list, num_genes, rr)

    return output_id_set, rr

def main():
    """ dump stableIDs for expressing genes in a study based on commonality """
    script_info = set_environment()
    option_parser, opts, args =\
            parse_command_line_parameters(**script_info)

    rr = RunRecord()

    if opts.sample1 is None:
        raise RuntimeError('No samples given')

    # These will hold the ids we want to intersect, etc
    sample1_ids = set()
    sample2_ids = set()
    sample3_ids = set()

    # Get all the genes and build ensembl ID sets
    session = _create_session()
    sample1_genes, rr = getExpressedGenes(session, opts.sample1,
            opts.sample1_type, opts.m1, opts.sample_extremes, rr)
    for gene in sample1_genes:
        sample1_ids.add(gene.ensembl_id)
    session.close()

    session = _create_session()
    sample2_genes, rr = getExpressedGenes(session, opts.sample2,
            opts.sample2_type, opts.m2, opts.sample_extremes, rr)
    for gene in sample2_genes:
        sample2_ids.add(gene.ensembl_id)
    session.close()

    if opts.sample3 is not None:
        session = _create_session()
        sample3_genes, rr = getExpressedGenes(session, opts.sample3,
                opts.sample3_type, opts.m3, opts.sample_extremes, rr)
        for gene in sample3_genes:
            sample3_ids.add(gene.ensembl_id)
        session.close()

    # Find the IDs we're interested in based on sample relationship
    comparison_type = opts.comparison_type.split(':')[0]
    output_id_set, rr = compare_ids(sample1_ids, sample2_ids, sample3_ids,
            comparison_type, rr)

    # Narrow search if needed by top genes
    if opts.num_genes is not None:
        output_id_set, rr = restrict_by_num_genes(output_id_set, opts.num_genes,
                opts.expression_sample1, opts.expression_sample2,
                opts.expression_sample3, opts.favoured_expression_sample, rr)


    rr.addInfo('gene_overlap', 'Total genes in output', len(output_id_set))
    # now save to file
    outfile = open(opts.genefile, 'w')
    # Add the standard header so that we can import with add_expression_db.py
    outfile.write('gene\n')
    for id in output_id_set:
        outfile.write(str(id)+'\n')
    outfile.close()
    rr.display()

if __name__ == '__main__':
    main()
