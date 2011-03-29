#!/usr/bin/env python
import os, re
from cogent.parse.psl import MinimalPslParser
from cogent.parse.fastq import MinimalFastqParser
from chippy.parse.light_seq import LightSeq

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"
__version__ = '0.1'

file_end = re.compile(r'\.fastq$')

def get_corrupt_seq_names(psl_name, test_run):
    psl_parser = MinimalPslParser(psl_name)
    psl_parser.next()
    psl_parser.next()
    
    num = 0
    contaminated_info = {}
    for record in psl_parser:
        # get the contaminated sequence name and the index where contamination
        # starts and store it as a key-->value pair
        contaminated_info[record[9]] = record[11]
        num += 1
        if test_run and num >= 1000:
            break
    return contaminated_info

def write_pristine(fastq_name, not_pristine, run_record, test_run):
    num = 0
    num_too_short = 0
    num_pristine = 0
    num_contaminated = 0
    pristine_name = file_end.sub('_pristine.fastq', fastq_name)
    contaminated_name = file_end.sub('_contaminated.fastq', fastq_name)
    
    if test_run:
        print pristine_name
        print contaminated_name
    else:
        outfile_pristine = open(pristine_name, 'w')
        outfile_contaminated = open(contaminated_name, 'w')
    
    seq_object = LightSeq()
    for name, seq, qual in MinimalFastqParser(open(fastq_name)):
        num += 1

        seq_object(seq, name, qual)

        if test_run:
            print fastq_formatted
            if num > 100:
                break

        try:
            start = not_pristine[name]

            # since we trim everything after the adapter sequence start,
            # there is little value in keeping a sequence which is going to be
            # less than 35 bp long.
            if start < 35:
                num_too_short += 1
                continue

            num_contaminated += 1
            seq_object = seq_object[:start]
            fastq_formatted = seq_object.toFastq()
            outfile_contaminated.write(fastq_formatted + '\n')

        # if the sequence doesnt exist in output from blat, it must not be
        # contaminated with the adapter
        except KeyError:
            num_pristine += 1
            fastq_formatted = seq_object.toFastq()
            outfile_pristine.write(fastq_formatted + '\n')

    if not test_run:
        outfile_contaminated.close()
        outfile_pristine.close()
    
    run_record.addMessage(program_name='pristine_seqs',
                error_type='stdout', message='Pristine Seqs',
                value=num)
    run_record.addMessage(program_name='pristine_seqs',
                error_type='stdout',
                message='Sequences were contaminated with adapter, but still kept',
                value=num_contaminated)
    run_record.addMessage(program_name='pristine_seqs',
                error_type='stdout',
                message='Sequences were discarded as too small',
                value=num_too_short)
    return run_record

def main(input_psl_file, input_file, run_record, test_run):
    """identifies reads not to be written, then writes everything else"""
    not_pristine = get_corrupt_seq_names(input_psl_file, test_run)
    run_record = write_pristine(input_file, not_pristine, run_record, test_run)
    return run_record

if __name__ == "__main__":
    from cogent.util.misc import parse_command_line_parameters
    from optparse import make_option

    script_info = {}
    descr = "Write seqs without primer/adapter match to separate file."
    script_info['brief_description']= descr
    script_info['script_description'] = descr
    script_info['version'] = '0.1.alpha'
    script_info['script_usage']=[]
    script_info['script_usage'].append(
        ("Example 1","""Test run of write pristine:""",
        """python pristine_seqs.py -p <somefile.psl> -i <seqs> -o <outfile_root> -t"""))
    
    script_info['help_on_no_arguments'] = True
    script_info['required_options'] = [
        make_option('-p','--input_psl_file',
                    help='The input psl file from blat'),
        make_option('-i','--input_file',
                    help='The input fastq sequence file'),
        ]
    
    script_info['optional_options'] = [\
        make_option('-t','--test_run', action='store_true',
                    dest='test_run', default = False,
                    help='Dry run without writing any data'
                    +'[default: %default]'),
                    ]
    
    parser, opts, args = parse_command_line_parameters(**script_info)
    
    main(opts.input_psl_file, opts.input_file, opts.test_run)
