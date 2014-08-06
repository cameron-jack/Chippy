#!/usr/bin/env python
"""
    Reads a WIG file of mRNA expression and converts it to
    ChipPy's expression file format (exp)
"""
from __future__ import division

import os, sys, warnings, gzip
warnings.filterwarnings('ignore', 'Not using MPI as mpi4py not found')
sys.path.extend(['..', '../src'])

from cogent.util.progress_display import display_wrap

from chippy.express import db_query
from chippy.util.run_record import RunRecord
from chippy.util.command_args import Args
from chippy.util.util import run_command

import numpy

__author__ = 'Cameron Jack'
__copyright__ = 'Copyright 2011-2013, Gavin Huttley, Cameron Jack, Anuj Pahwa'
__credits__ = ['Cameron Jack']
__license__ = 'GPL'
__maintainer__ = 'Cameron Jack'
__email__ = 'cameron.jack@anu.edu.au'
__status__ = 'Pre-release'
__version__ = '0.2'

script_info = {}
script_info['title'] = 'Convert WIG to EXP file'
script_info['script_description'] = 'Reads a WIG file of mRNA expression and '+\
        "converts it to ChipPy's expression file format. Handles gzipped files."
script_info['brief_description'] = 'Convert WIG to flat file'
script_info['version'] = __version__
script_info['authors'] = __author__
script_info['output_description']= 'File with gene and expr values'
pos_args = ['db_path']
req_args = ['wig']
opt_args = ['max_chrom_size', 'chr_prefix', 'exp']

script_info['args'] = Args(required_args=req_args, optional_args=opt_args,
    positional_args=pos_args)
script_info['required_options'] = script_info['args'].getReqCogentOpts()
script_info['optional_options'] = script_info['args'].getOptCogentOpts()

def get_gene_scores_from_chrom(chrom_array, chrom, all_genes, genes_by_chrom,
        genes_scores):
    """ for genes in the current chrom, slice chrom_array and sum the scores """
    try:
        id_list = genes_by_chrom[chrom]
    except KeyError:
        id_list = []
    for id in id_list:
        gene = all_genes[id]
        score = sum(chrom_array[gene.start:gene.end])
        genes_scores[gene.ensembl_id] = score

@display_wrap
def main(ui=None):
    """
        1) Get all protein coding genes from DB.
        2) Read WIG file and if a count record is in a gene then add
            to its total
        3) Write out genes and expr values
    """
    rr = RunRecord('expr_wig_to_exp')
    rr.addCommands(sys.argv)

    args = script_info['args'].parse(window_title='Expression WIG to EXP')
    chrom_size = args.max_chrom_size
    prefix = args.chr_prefix

    session = db_query.make_session(args.db_path)
    genes = db_query.get_gene_entries(session)

    all_genes = {} # genes indexed by ensembl_id
    genes_by_chrom = {} # chrom: list(gene_id)
    genes_scores = {} # each gene has an expression score
    for gene in genes:
        if not gene.chrom in genes_by_chrom.keys():
            genes_by_chrom[gene.chrom] = []
        genes_by_chrom[gene.chrom].append(gene.ensembl_id)
        genes_scores[gene.ensembl_id] = 0
        all_genes[gene.ensembl_id] = gene

    wig_fn = args.wig
    if wig_fn.endswith('.gz'):
        wig_file = gzip.GzipFile(wig_fn, 'rb')
    else:
        try:
            wig_file = open(wig_fn, 'r')
        except IOError:
            rr.dieOnCritical('Could not open file', wig_fn)

    # get total lines in wig for pacing the progress bar
    if not wig_fn.endswith('.gz'):
        command = 'wc -l ' + wig_fn
        returncode, stdout, stderr = run_command(command)
        if returncode:
            rr.addWarning('could not run wc to count WIG lines', 'error')
            total_lines = 1
        else:
            total_lines = int(stdout.strip().split(' ')[0])
            rr.addInfo('total lines in '+wig_fn, total_lines)

    # Read each piece of the file into an artificial chromosome (Numpy array)
    # and slice out the gene regions that we have for each gene in that chrom

    chrom_array = numpy.zeros(chrom_size, dtype=numpy.float32)

    current_chrom = None
    for i, line in enumerate(wig_file):
        if i % 100 == 0:
            msg = 'Reading wiggle entries [' + str(i) +\
                  ' / ' + str(total_lines) + ']'
            progress = (float(i)/float(total_lines))
            ui.display(msg=msg, progress=progress)

        if line.startswith('track'):
            continue
        elif line.startswith('fixed'):
            # fixedStep chrom=chr10 start=56001 step=20 span=20
            step_type = 'fixed'
            step_parts = line.split(' ')
            step = [val.strip('step=').strip() \
                    for val in step_parts if val.startswith('step')][0]
            span = [val.strip('span=').strip() \
                    for val in step_parts if val.startswith('span')][0]
            chrom = [val.strip('chrom='+prefix).strip() \
                     for val in step_parts if val.startswith('chrom')][0]

            if chrom == 'M':
                chrom = 'MT'

            if current_chrom is None:
                current_chrom = chrom
            elif current_chrom != chrom: # Empty chrom_array into genes
                get_gene_scores_from_chrom(chrom_array, chrom, all_genes,
                        genes_by_chrom, genes_scores)
                current_chrom = chrom
                chrom_array[:] = 0

            start = [val.strip('start=').strip() \
                     for val in step_parts if val.startswith('start')][0]
            pos = int(start)
            step = int(step)
            span = int(span)
        elif line.startswith('variable'):
            step_type = 'variable'
            step_parts = line.split(' ')
            chrom = [val.strip('chrom='+prefix).strip() \
                    for val in step_parts if val.startswith('chrom')][0]

            if chrom == 'M':
                chrom = 'MT'

            if current_chrom is None:
                current_chrom = chrom
            elif current_chrom != chrom: # Empty chrom_array into genes
                get_gene_scores_from_chrom(chrom_array, chrom, all_genes,
                        genes_by_chrom, genes_scores)
                current_chrom = chrom
                chrom_array[:] = 0
        else:
            if step_type == 'fixed':
                chrom_array[pos] = float(line.strip())
                pos += step
            else: #step_type == 'variable'
                if '\t' in line:
                    line_parts = line.split('\t')
                else:
                    line_parts = line.split(' ')
                chrom_array[int(line_parts[0])] = float(line_parts[1].strip())

    # empty chrom_array into genes_score from the final section
    get_gene_scores_from_chrom(chrom_array, chrom, all_genes,
            genes_by_chrom, genes_scores)

    # output genes and scores
    if args.exp:
        out_fn = args.exp
    else:
        if '.gz' in wig_fn:
            wig_fn = '.'.join(wig_fn.split('.')[:-1])
        out_fn = '.'.join(wig_fn.split('.')[:-1]) # cut off wig extension
        out_fn += '.exp' # add .exp extension

    with open(out_fn, 'w') as out:
        out.write('gene\texp\n') # header
        for id in genes_scores.keys():
            out.write(id + '\t' + str(genes_scores[id]) + '\n')
        out.close()

if __name__ == '__main__':
    main()