#!/usr/bin/env python
"""dumps a RegionCollection compressed file"""
from __future__ import division

import os, sys, glob, warnings
warnings.filterwarnings('ignore', 'Not using MPI as mpi4py not found')
sys.path.extend(['..', '../src'])

from chippy.core.count_tags import counts_for_genes
from chippy.express import db_query
from chippy.express.util import sample_types
from chippy.util.run_record import RunRecord
from chippy.util.command_args import Args

__author__ = 'Cameron Jack, Gavin Huttley'
__copyright__ = 'Copyright 2011-2013, Gavin Huttley, Cameron Jack, Anuj Pahwa'
__credits__ = ['Gavin Huttley', 'Cameron Jack']
__license__ = 'GPL'
__maintainer__ = 'Cameron Jack'
__email__ = 'cameron.jack@anu.edu.au'
__status__ = 'Pre-release'
__version__ = '0.1'

script_info = {}
script_info['title'] = 'Saves feature counts'
script_info['script_description'] = 'Saves counts data for TSS, ' +\
        'Exon-Intron, Intron-Exon or Gene-3prime boundaries ' +\
        'for a given window. Reads count info from an indexed BAM ' +\
        'or a BED file'
script_info['brief_description'] = 'Extracts counts from sequenced regions'
script_info['version'] = __version__
script_info['authors'] = __author__
script_info['output_description']= 'Generates a Pickle file or a ' +\
        'gzipped tab-delimited file that can be used for plotting ' +\
        'of subsets of genes. Can also output to BEDgraph.'
pos_args = ['db_path']
req_args = ['sample', 'sample_type', 'feature_type',
        'BAMorBED',  'collection']
opt_args = ['overwrite', 'tab_delimited', 'max_read_length', 'chr_prefix',
        'count_max_length', 'window_upstream', 'window_downstream',
        'multitest_signif_val', 'include_target',
        'exclude_target', 'make_bedgraph']

script_info['args'] = Args(required_args=req_args, optional_args=opt_args,
    positional_args=pos_args)
script_info['required_options'] = script_info['args'].req_cogent_opts
script_info['optional_options'] = script_info['args'].opt_cogent_opts

def get_collection(session, sample_name, sample_type, feature_type, BAMorBED,
        chr_prefix, window_upstream, window_downstream,
        multitest_signif_val, collection_fn, overwrite,
        tab_delimited, include_target=None, exclude_target=None,
        bedgraph=False):
    """
        builds and writes a collection of counts and expression for
        feature_type in given sample genes.
    """
    rr = RunRecord('get_collection')

    if not os.path.exists(collection_fn) or overwrite:
        bedgraph_fn = None
        if bedgraph:
            bedgraph_fn = collection_fn.split('.')[0] + '.bedgraph'

        if sample_type == sample_types['exp_absolute'] or \
                sample_type == sample_types['exp_diff']:
            data_collection = counts_for_genes(session,
                    sample_name, sample_type, feature_type, BAMorBED,
                    chr_prefix, window_upstream, window_downstream,
                    include_target, exclude_target, bedgraph_fn,
                    multitest_signif_val=multitest_signif_val)
        else:
            rr.dieOnCritical('Experiment type not supported', sample_type)

        if data_collection is not None:
            data_collection.writeToFile(collection_fn, as_table=tab_delimited,
                    compress_file=True)
        else:
            rr.dieOnCritical('No data collection was returned', 'Failed')
    else:
        print 'Existing output at', collection_fn

def main():
    """
        Returns a pickle of size window_start to window_finish containing
        chromatin mapping averages per base, one per gene, ranked by
        expression.
    """
    rr = RunRecord('export_counts')
    rr.addCommands(sys.argv)

    args = script_info['args'].parse()
    if args.sample is None:
        rr.dieOnCritical('No samples provided', 'Failed')

    db_name = str(args.db_path).split('/')[-1]

    sample_name = args.sample.split(' : ')[0]
    print 'Loading counts data for', sample_name
    sample_type = sample_types[args.sample_type]

    include_name = None
    exclude_name = None
    if args.include_target:
        include_name = args.include_target.split(' : ')[0]
        rr.addInfo('include gene targets', include_name)

    if args.exclude_target:
        exclude_name = args.exclude_target.split(' : ')[0]
        rr.addInfo('exclude gene targets', exclude_name)

    if (args.multitest_signif_val is not None) and not \
            (-1 <= args.multitest_signif_val <= 1):
        rr.dieOnCritical('Multitest_signif_val should be -1, 0, 1',
                args.multitest_signif_val)

    session = db_query.make_session(args.db_path)

    if args.chr_prefix != '':
        # If it writes nothing then we'll fail to load the table again
        rr.addInfo('BAM/BED chromosome prefix given', args.chr_prefix)

    window_upstream = args.window_upstream
    assert window_upstream > 0, \
            'upstream window must be of at least size 1 bp'
    window_downstream = args.window_downstream
    assert window_downstream > 0, \
            'downstream window must be of at least size 1 bp'

    if sample_type in [sample_types['exp_absolute'],\
            sample_types['exp_diff']]:
        get_collection(session, sample_name, sample_type, args.feature_type,
                args.BAMorBED, args.chr_prefix, window_upstream,
                window_downstream, args.multitest_signif_val,
                args.collection, args.overwrite,
                args.tab_delimited, include_name,
                exclude_name, bedgraph=args.make_bedgraph)
    else:
        print 'Other options not defined yet, choose from', \
                sample_types['exp_absolute'], 'or', sample_types['exp_diff']

    session.close()
    rr.display()

if __name__ == '__main__':
    main()

