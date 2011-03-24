#!/usr/bin/env python
import time
from os.path import basename, dirname, join

from cogent.parse.fastq import MinimalFastqParser
from cogent.util.misc import parse_command_line_parameters
from optparse import make_option

from chippy.util.util import create_path
from chippy.parse.fastq import FastqParser

def run(input_file, save_dir, output_file, minimum_length, rewrite, run_record, test_run):
    if not test_run:
        assert '/' not in output_file, \
                'Do not put directory path in output_file'
        
        create_path(save_dir)
        outfile_fasta = open(join(save_dir, output_file), 'w')
        if rewrite:
            outfile_fastq_fn = basename(output_file).split('.')[0] + '_trimmed.fastq'
            outfile_fastq = open(join(save_dir, outfile_fastq_fn), 'w')
    
    i = 0
    num_too_small = 0
    for seq in FastqParser(input_file):
        i += 1
        if len(seq) < minimum_length:
            num_too_small += 1
            continue

        data_fasta = seq.toFasta() + '\n'
        if rewrite:
            data_fastq = seq.toFastq() + '\n'
        
        if not test_run:
            outfile_fasta.write(data_fasta)
            if rewrite:
                outfile_fastq.write(data_fastq)
        else:
            print data_fasta

        if test_run and i > 10000:
            break

    if not test_run:
        outfile_fasta.close()
        if rewrite:
            outfile_fastq.close()
    
    run_record.addMessage(program_name='fastq_to_fasta',
                error_type='stdout', message='Sequences read', value=i)
    run_record.addMessage(program_name='fastq_to_fasta',
            error_type='stdout', message='Sequences discarded as too small',
            value=num_too_small)
    
    return run_record

# 
script_info = {}
descr = "Convert fastq file to a fasta file after quality checks on fastq file and "\
        "only choosing sequences that meet specified criteria"
script_info['brief_description']= descr
script_info['script_description'] = descr
script_info['version'] = '0.1.alpha'
script_info['script_usage']=[]
script_info['script_usage'].append(
    ("Example 1","""General Usage:""",
    """python fastq_to_fasta.py -i <seqs.fastq> -o <seqs.fasta>"""))

script_info['help_on_no_arguments'] = True
script_info['required_options'] = [
    make_option('-i','--input_file',
                help='The input fastq sequence file'),
    make_option('-o','--output_file',
                help='The file to write fasta seqs to'),
    make_option('-S', '--save_dir',
                help='path to save all files')
    ]

script_info['optional_options'] = [\
    make_option('-t','--test_run', action='store_true',
                dest='test_run', default = False,
                help='Dry run without writing any data'
                +'[default: %default]'),
    make_option('-l', '--minimum_length', type='int',
                default=35,
                help='minimum length of sequences to write [default: %default]'),
    make_option('-q', '--rewrite_fastq', action='store_true',
                help='Rewrite input fastq file after removing discarded sequences')
                ]

def main():
    from fastq_to_mapped import RunRecord
    run_record = RunRecord()
    start = time.time()
    parser, opts, args = parse_command_line_parameters(**script_info)
    run_record = run(opts.input_file, opts.save_dir, opts.output_file,
        opts.minimum_length, opts.rewrite_fastq, run_record, opts.test_run)
    end = time.time()
    run_record.addMessage(program_name='fastq_to_fasta',
            error_type='stdout', message='Time taken (mins)',
            value= ((end-start)/60.))
    return run_record

if __name__ == "__main__":
    from cogent import LoadTable
    run_record = main()
    run_record.display()

