#!/usr/bin/env python
"""dumps a RegionCollection compressed file"""
from __future__ import division

import os, sys, warnings
warnings.filterwarnings('ignore', 'Not using MPI as mpi4py not found')
sys.path.extend(['..', '../src'])

from chippy.core.count_tags import counts_for_genes
from chippy.express import db_query
from chippy.util.run_record import RunRecord
from chippy.util.command_args import Args

__author__ = 'Cameron Jack, Gavin Huttley'
__copyright__ = 'Copyright 2011-2013, Gavin Huttley, Cameron Jack, Anuj Pahwa'
__credits__ = ['Gavin Huttley', 'Cameron Jack']
__license__ = 'GPL'
__maintainer__ = 'Cameron Jack'
__email__ = 'cameron.jack@anu.edu.au'
__status__ = 'Pre-release'
__version__ = '0.2'

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
req_args = ['expr_sample', 'feature_type', 'BAMorBED',  'collection']
opt_args = ['overwrite', 'tab_delimited', 'max_read_length', 'chr_prefix',
        'count_max_length', 'window_upstream', 'window_downstream',
        'multitest_signif_val', 'include_target',
        'exclude_target', 'make_bedgraph']

script_info['args'] = Args(required_args=req_args, optional_args=opt_args,
    positional_args=pos_args)
script_info['required_options'] = script_info['args'].req_cogent_opts
script_info['optional_options'] = script_info['args'].opt_cogent_opts

def get_collection(session, sample_name, feature_type, BAMorBED,
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

        data_collection = counts_for_genes(session, sample_name, feature_type,
                BAMorBED, chr_prefix, window_upstream, window_downstream,
                include_target, exclude_target, bedgraph_fn,
                multitest_signif_val=multitest_signif_val)

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

    args = script_info['args'].parse(use_scrollbars=True,
            use_save_load_button=True, window_title='Export Counts')

    session = db_query.make_session(args.db_path)

    sample_name = args.expr_sample
    print 'Loading counts data for', sample_name

    include_name = None
    exclude_name = None
    if args.include_target:
        include_name = args.include_target
        rr.addInfo('include gene targets', include_name)

    if args.exclude_target:
        exclude_name = args.exclude_target
        rr.addInfo('exclude gene targets', exclude_name)

    if (args.multitest_signif_val is not None) and not \
            (-1 <= args.multitest_signif_val <= 1):
        rr.dieOnCritical('Multitest_signif_val should be -1, 0, 1',
                args.multitest_signif_val)

    if args.chr_prefix != '':
        # If it writes nothing then cogent.Table fails because it's fragile
        rr.addInfo('BAM/BED chromosome prefix given', args.chr_prefix)

    window_upstream = args.window_upstream
    assert window_upstream > 0, \
            'upstream window must be of at least size 1 bp'
    window_downstream = args.window_downstream
    assert window_downstream > 0, \
            'downstream window must be of at least size 1 bp'

    get_collection(session, sample_name, args.feature_type, args.BAMorBED,
            args.chr_prefix, window_upstream, window_downstream,
            args.multitest_signif_val, args.collection, args.overwrite,
            args.tab_delimited, include_name, exclude_name,
            bedgraph=args.make_bedgraph)

    session.close()
    rr.display()

if __name__ == '__main__':
    main()

