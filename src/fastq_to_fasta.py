#!/usr/bin/env python

from os.path import basename, dirname, join
from cogent.parse.fastq import MinimalFastqParser
from util import create_path

def make_seq(seq, name, qual):
    """over-ride default, slower, sequence constructor"""
    return name, seq

def run(input_file, output_file, minimum_length, rewrite, test_run):
    if not test_run:
        outfile_fasta = open(output_file, 'w')
        outfile_directory = dirname(output_file)
        create_path(outfile_directory)
        if rewrite:
            outfile_fastq_fn = basename(output_file).split('.')[0] + '_trimmed.fastq'
            outfile_fastq = open(join(outfile_directory,outfile_fastq_fn), 'w')
    
    i = 0
    num_too_small = 0
    #print '\nShort Sequences:\n'
    for seq in MinimalFastqParser(input_file):
        i += 1
        if len(seq) < minimum_length:
            #print '%s:\t%d\n' % (seq.Name,len(seq))
            num_too_small += 1
            continue

        #if i % 100000 == 0:
            #print seq.Name
            #print i

        data_fasta = seq.toFasta() + '\n'
        if rewrite:
            data_fastq = seq.toFastq() + '\n'
        if not test_run:
            outfile_fasta.write(data_fasta)
            if rewrite:
                outfile_fastq.write(data_fastq)
        else:
            print data

        if test_run and i > 10000:
            break

    if not test_run:
        outfile_fasta.close()
        if rewrite:
            outfile_fastq.close()

    print 'Success!!\n%d Sequences were read' % i
    print '%d Sequences were discarded as too small' % num_too_small

if __name__ == "__main__":
    from cogent.util.misc import parse_command_line_parameters
    from optparse import make_option

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
                    help='The file to write fasta seqs to')
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

    parser, opts, args = parse_command_line_parameters(**script_info)
    run(opts.input_file, opts.output_file, opts.minimum_length, opts.rewrite_fastq, opts.test_run)