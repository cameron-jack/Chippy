#!/usr/bin/env python
"""implements the full workflow of processing sequences
- converts to fasta in prep for using blat
- uses blat to find adapter sequences
- produces trimmed sequences without adpaters
- maps the clean sequences using bowtie
"""
import sys, re, os, shutil
sys.path.extend(['..', '../src'])

import numpy

from cogent.util.misc import parse_command_line_parameters
from optparse import make_option

from chippy.prep import reduce, pristine_seqs, command_line, fastq_to_fasta, \
     pristine_paired
from chippy.prep.mapped_files import MappedFiles
from chippy.util.run_record import RunRecord
from chippy.util.util import create_path


script_info = {}
descr = "Draft snp discovery process"
script_info['brief_description']= descr
script_info['script_description'] = descr
script_info['version'] = '1e-3.alpha'
script_info['script_usage']=[]
# script_info['script_usage'].append(
#     ("Example 1","""General Usage:""",
#     """python fastq_to_fasta.py -i <seqs.fastq> -o <seqs.fasta>"""))

script_info['help_on_no_arguments'] = True
script_info['required_options'] = [
    make_option('-1','--input_file_1',
                help='The 1st input fastq sequence file'),
    make_option('-2','--input_file_2',
                help='The 2nd input fastq sequence file'),
    make_option('-S', '--save_dir', help='path to save all files'),
    make_option('--adapters', help='path to the Illumina adapters'),
    make_option('--index', help='path to the bwa aligner index')
]

script_info['optional_options'] = [\
    make_option('-n', '--num_threads', type='int',
                default = 6,
                help='Number of threads to use [default: %default]'),
    make_option('-m', '--mem_usage', type='int',
                default = 7000000000,
                help='memory usage for samtools sort [default: %default]'),
    make_option('-t', '--test_run', action='store_true',
                dest= 'test_run', default = False,
                help='Dry run without writing any data'
                +' [default: %default]'),
    make_option('-s', '--sample_name', type='string', default = '',
                help='specify sample name in annotated bam header'
                +' [default: %default]'),
    make_option('-w', '--work_dir', type='string', default = '',
                help='specify temporary working directory'
                +' [default: %default]'),
    make_option('-b', '--begin', type='int', default = 1,
                help='begin at stage # [default: %default]'),
    make_option('-e', '--end', type='int', default = 7,
                help='end at stage # [default: %default]'),
    make_option('-D', '--delete', action='store_true', default = False,
                help='Deletes the working dir at the end of the run'
                +' [default: %default]'),
    make_option('-r', '--reduce', action='store_true', default = False,
                help='Finish with ChIP-Seq reduction step'
                +' [default: %default]'),
    make_option('-I', '--Illumina_version', default=1.7,
                help='Illumina pipeline version number'
                +' [default: %default]'),
    make_option('--no_qual', action='store_true', default=False,
                help='Do not use any quality control'
                +' [default: %default]')
]

def main():
    option_parser, opts, args =\
        parse_command_line_parameters(**script_info)

    # make the assorted filenames

    rr = RunRecord() # the run audit tracker

    mapped_files_1 = MappedFiles(opts.input_file_1, opts.save_dir, opts.work_dir)
    filenames_1 = dict(fastq=mapped_files_1.adapterless_trimmed_fn,
                       pristine=mapped_files_1.pristine_fn, sai=mapped_files_1.sai_fn,
                       bam=mapped_files_1.bam_fn, sam=mapped_files_1.sam_fn)
    mapped_files_2 = MappedFiles(opts.input_file_2, opts.save_dir, opts.work_dir)
    filenames_2 = dict(fastq=mapped_files_2.adapterless_trimmed_fn,
                       pristine=mapped_files_2.pristine_fn, sai=mapped_files_2.sai_fn,
                       bam=mapped_files_2.bam_fn)

    if os.path.exists(filenames_1['bam']):
        rr.addInfo('fastq_to_mapped_bam', 'already exists, exiting', filenames_1['bam'])
        rr.display()
        sys.exit(0)

    if not opts.no_qual:
        if 1.3 >= opts.Illumina_version <= 1.7:

            ## do input file 1
            if opts.begin <= 1 and opts.end >= 1:
                rr = command_line.run_fastx_clip_and_trim(opts.adapters, opts.input_file_1,
                                                  filenames_1['fastq'], rr, opts.num_threads, opts.test_run)

            ## do input file 2
            if opts.begin <= 2 and opts.end >= 2:
                rr = command_line.run_fastx_clip_and_trim(opts.adapters, opts.input_file_2,
                                                  filenames_2['fastq'], rr, opts.num_threads, opts.test_run)

            ## produce the paired pristine seq files with matched reads
            if opts.begin <= 3 and opts.end >= 3:
                # rr = command_line.run_pristine_paired(
                rr = command_line.run_pristine_seesaw(opts.input_file_1,
                                              filenames_1['fastq'], filenames_2['fastq'],
                                              filenames_1['pristine'], filenames_2['pristine'],
                                              rr, opts.num_threads, opts.test_run)

        elif opts.Illumina_version >= 1.8: # Illumina pipeline 1.8 or higher
            # Don't worry about 5' adapter clipping
            # Though we might want to think about 3' adapters in future

            # Quality check with Sickle (UC Davis), discard single-ends residual
            quality_window = 20;
            min_length = 19;
            rr = command_line.run_sickle_pe(opts.input_file_1, opts.input_file_2,
                            filenames_1['pristine'], filenames_2['pristine'],
                            quality_window, min_length, rr, opts.test_run)

        else:
            rr.addError('fastq_to_mapped_bam', 'Invalid Illumina pipeline given. v1.3+ supported',
                        opts.Illumina_version)
            rr.display()
            sys.exit(0)

        ## now align each individually
        if opts.begin <= 4 and opts.end >= 4:
            rr = command_line.run_bwa_aln(opts.index,
                                      filenames_1['pristine'],
                                      filenames_1['sai'],
                                      opts.Illumina_version,
                                      rr, opts.num_threads, opts.test_run)

        if opts.begin <= 5 and opts.end >= 5:
            rr = command_line.run_bwa_aln(opts.index,
                                      filenames_2['pristine'],
                                      filenames_2['sai'],
                                      opts.Illumina_version,
                                      rr, opts.num_threads, opts.test_run)

        ## direct bwa_sampe to final sorted bam, bypassing sam
        if opts.begin <= 6 and opts.end >= 6:
            # rr = command_line.bwa_sampe_to_sorted_bam(opts.index,
            rr = command_line.bwa_sampe_to_local_bam_sort_index_move(opts.index,
                                      filenames_1['sai'], filenames_2['sai'],
                                      filenames_1['pristine'], filenames_2['pristine'],
                                      filenames_1['bam'],
                                      rr, opts.sample_name, opts.mem_usage, opts.test_run)

    else:
        ## read direct from fastq
        if opts.begin <= 4 and opts.end >= 4:
            rr = command_line.run_bwa_aln(opts.index,
                                      opts.input_file_1,
                                      filenames_1['sai'],
                                      opts.Illumina_version,
                                      rr, opts.num_threads, opts.test_run)

        if opts.begin <= 5 and opts.end >= 5:
            rr = command_line.run_bwa_aln(opts.index,
                                      opts.input_file_2,
                                      filenames_2['sai'],
                                      opts.Illumina_version,
                                      rr, opts.num_threads, opts.test_run)

        ## direct bwa_sampe to final sorted bam, bypassing sam
        if opts.begin <= 6 and opts.end >= 6:
            # rr = command_line.bwa_sampe_to_sorted_bam(opts.index,
            rr = command_line.bwa_sampe_to_local_bam_sort_index_move(opts.index,
                                      filenames_1['sai'], filenames_2['sai'],
                                      opts.input_file_1, opts.input_file_2,
                                      filenames_1['bam'],
                                      rr, opts.sample_name, opts.mem_usage, opts.test_run)

    if opts.reduce:
        # convert to SAM and reduce
        rr = command_line.convert_bam_to_sam(filenames_1['bam'],
                                             filenames_1['sam'],rr, opts.test_run)

        rr = reduce.run(infile_name=filenames_1['sam'],
        outdir=opts.save_dir, chroms='Do All', pval_cutoff=opts.pval_cutoff,
        limit=numpy.inf, run_record=rr, dry_run=opts.test_run)

    if opts.delete:
        os.remove(filenames_1['sai'])
        os.remove(filenames_2['sai'])
        os.remove(filenames_1['pristine'])
        os.remove(filenames_2['pristine'])
        os.remove(filenames_1['sam']) # Keep BAM

        ## index final bam
#    if opts.begin <= 7 and opts.end >= 7:
#        rr = command_line.index_bam(filenames_1['bam'], rr, opts.test_run)

    ## output audit
    rr.display()
    table = rr.getMessageTable()
    table.writeToFile(mapped_files_1.run_record_fn, sep='\t')

if __name__ == "__main__":
    main()