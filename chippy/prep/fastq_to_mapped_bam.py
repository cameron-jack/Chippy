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
    make_option('--index', help='path to the bwa aligner index'),
]

script_info['optional_options'] = [\
    make_option('-n','--num_threads', type='int',
                default = 6,
                help='Number of threads to use [default: %default]'),
    make_option('-m','--mem_usage', type='int',
                default = 7000000000,
                help='memory usage for samtools sort [default: %default]'),
    make_option('-t','--test_run', action='store_true',
                dest='test_run', default = False,
                help='Dry run without writing any data'
                +'[default: %default]'),
    make_option('-s','--sample_name', type='string', default = '',
                help='specify sample name in annotated bam header'
                +'[default: %default]'),
    make_option('-w','--work_dir', type='string', default = '',
                help='specify temporary working directory'
                +'[default: %default]'),
    make_option('-b','--begin', type='int', default = 1,
                help='begin at stage # [default: %default]'),
    make_option('-e','--end', type='int', default = 7,
                help='end at stage # [default: %default]'),
    make_option('-D','--delete', action='store_true', default = False,
                help='Deletes the working dir at the end of the run'
                +'[default: %default]'),
]

def main():
    option_parser, opts, args =\
        parse_command_line_parameters(**script_info)

    # make the assorted filenames

    rr = RunRecord() # the run audit tracker

    mapped_files_1 = MappedFiles(opts.input_file_1, opts.save_dir, opts.work_dir)
    filenames_1 = dict(fastq=mapped_files_1.adapterless_trimmed_fn,
                       pristine=mapped_files_1.pristine_fn, sai=mapped_files_1.sai_fn,
                       bam=mapped_files_1.bam_fn)
    mapped_files_2 = MappedFiles(opts.input_file_2, opts.save_dir, opts.work_dir)
    filenames_2 = dict(fastq=mapped_files_2.adapterless_trimmed_fn,
                       pristine=mapped_files_2.pristine_fn, sai=mapped_files_2.sai_fn,
                       bam=mapped_files_2.bam_fn)

    if os.path.exists(filenames_1['bam']):
        rr.addInfo('fastq_to_mapped_bam', 'already exists, exiting', filenames_1['bam'])
        rr.display()
        exit(0)

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
#        rr = command_line.run_pristine_paired(
        rr = command_line.run_pristine_seesaw(opts.input_file_1,
                                              filenames_1['fastq'], filenames_2['fastq'],
                                              filenames_1['pristine'], filenames_2['pristine'],
                                              rr, opts.num_threads, opts.test_run)

    if opts.delete:
        os.remove(filenames_1['fastq'])
        os.remove(filenames_2['fastq'])

    ## now align each individually
    if opts.begin <= 4 and opts.end >= 4:
        rr = command_line.run_bwa_aln(opts.index,
                                      filenames_1['pristine'],
                                      filenames_1['sai'], rr, opts.num_threads, opts.test_run)

    if opts.begin <= 5 and opts.end >= 5:
        rr = command_line.run_bwa_aln(opts.index,
                                      filenames_2['pristine'],
                                      filenames_2['sai'], rr, opts.num_threads, opts.test_run)

    ## direct bwa_sampe to final sorted bam, bypassing sam
    if opts.begin <= 6 and opts.end >= 6:
#        rr = command_line.bwa_sampe_to_sorted_bam(opts.index,
        rr = command_line.bwa_sampe_to_local_bam_sort_index_move(opts.index,
                                                                 filenames_1['sai'], filenames_2['sai'],
                                                                 filenames_1['pristine'], filenames_2['pristine'],
                                                                 filenames_1['bam'],
                                                                 rr, opts.sample_name, opts.mem_usage, opts.test_run)

    if opts.delete:
        os.remove(filenames_1['sai'])
        os.remove(filenames_2['sai'])
        os.remove(filenames_1['pristine'])
        os.remove(filenames_2['pristine'])

        ## index final bam
#    if opts.begin <= 7 and opts.end >= 7:
#        rr = command_line.index_bam(filenames_1['bam'], rr, opts.test_run)

    ## output audit
    rr.display()
    table = rr.getMessageTable()
    table.writeToFile(mapped_files_1.run_record_fn, sep='\t')

if __name__ == "__main__":
    main()