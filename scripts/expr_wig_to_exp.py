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
opt_args = []

script_info['args'] = Args(required_args=req_args, optional_args=opt_args,
    positional_args=pos_args)
script_info['required_options'] = script_info['args'].getReqCogentOpts()
script_info['optional_options'] = script_info['args'].getOptCogentOpts()

def genes_in_block(pos_block, genes_by_chrom, chrom):
    """ get ids of any genes involved in current block """
    gene_ids = []
    if chrom not in genes_by_chrom.keys():
        return gene_ids
    for gene in genes_by_chrom[chrom]:
        if gene.start >= pos_block[0]:
            if gene.end <= pos_block[-1]:
                gene_ids.append(gene.ensembl_id)
    return gene_ids

def add_scores_genes(score_block, pos_block, ids, all_genes, genes_scores):
    """ for identified genes involved in a region add to their scores """
    for id in ids:
        gene = all_genes[id]
        start = gene.start
        end = gene.end

        scores =  numpy.array(score_block)
        for i, pos in enumerate(pos_block):
            if start <= pos <= end:
                genes_scores[id] += scores[i]

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

    session = db_query.make_session(args.db_path)
    genes = db_query.get_gene_entries(session)

    all_genes = {} # genes indexed by ensembl_id
    genes_by_chrom = {} # to save time in parsing positions
    genes_scores = {} # each gene has an expression score
    for gene in genes:
        if not gene.chrom in genes_by_chrom.keys():
            genes_by_chrom[gene.chrom] = []
        genes_by_chrom[gene.chrom].append(gene)
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

    # it is too slow to just check each position for its inclusion in a gene
    # We need to read 'blocks' of position data and then check these against
    # our gene coords
    pos_block = []
    score_block = []
    for i, line in enumerate(wig_file):
        if i % 100 == 0:
            msg = 'Reading wiggle entries [' + str(i) +\
                  ' / ' + str(total_lines) + ']'
            progress = (float(i)/float(total_lines))
            ui.display(msg=msg, progress=progress)

        if line.startswith('track'):
            continue
        elif line.startswith('fixed'):
            # empty pos and score blocks into genes_scores as appropriate
            if len(pos_block) > 0:
                ids = genes_in_block(pos_block, genes_by_chrom, chrom)
                if len(ids) > 0: # add scores to appropriate genes
                    add_scores_genes(score_block, pos_block, ids, all_genes,
                            genes_scores)
                pos_block = []
                score_block = []

            # fixedStep chrom=chr10 start=56001 step=20 span=20
            step_type = 'fixed'
            step_parts = line.split(' ')
            step = [val.strip('step=').strip() \
                    for val in step_parts if val.startswith('step')][0]
            span = [val.strip('span=').strip() \
                    for val in step_parts if val.startswith('span')][0]
            chrom = [val.strip('chrom=').strip() \
                     for val in step_parts if val.startswith('chrom')][0]
            start = [val.strip('start=').strip() \
                     for val in step_parts if val.startswith('start')][0]
            pos = int(start)
            step = int(step)
        elif line.startswith('variable'):
            # empty pos and score blocks into genes_scores as appropriate
            if len(pos_block) > 0:
                ids = genes_in_block(pos_block, genes_by_chrom, chrom)
                if len(ids) > 0: # add scores to appropriate genes
                    add_scores_genes(score_block, pos_block, ids, all_genes,
                        genes_scores)
                pos_block = []
                score_block = []

            step_type = 'variable'
            step_parts = line.split(' ')
            chrom = [val.strip('chrom=').strip() \
                     for val in step_parts if val.startswith('chrom')][0]
        else:
            if step_type == 'fixed':
                pos_block.append(pos)
                score_block.append(float(line.strip()))
                pos += step

            else: #step_type == 'variable'
                if '\t' in line:
                    line_parts = line.split('\t')
                else:
                    line_parts = line.split(' ')
                pos_block.append(int(line_parts[0]))
                score_block.append(float(line_parts[1].strip()))

    # empty pos and score blocks into genes_score from the final section
    if len(pos_block) > 0:
        ids = genes_in_block(pos_block, genes_by_chrom, chrom)
        if len(ids) > 0: # add scores to appropriate genes
            add_scores_genes(score_block, pos_block, ids, all_genes,
                    genes_scores)
            pos_block = []
            score_block = []

    # output genes and scores
    if '.gz' in wig_fn:
        wig_fn = ''.join(wig_fn.split('.')[:-1])
    out_fn = ''.join(wig_fn.split('.')[:-1]) # cut off wig extension
    out_fn += '.exp' # add .exp extension

    with open(out_fn, 'w') as out:
        out.write('gene\texp\n') # header
        for id in genes_scores.keys():
            out.write(id + '\t' + str(genes_scores[id]) + '\n')
        out.close()

if __name__ == '__main__':
    main()