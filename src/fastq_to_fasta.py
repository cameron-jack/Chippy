from parse_fastq import FastqParser

def make_seq(name, seq, qual):
    """over-ride default, slower, sequence constructor"""
    return name, seq

def run(input_file, output_file, minimum_length, test_run, verbose):
    if not test_run:
        outfile = open(output_file, 'w')
    i = 0
    num_too_small = 0
    for name, seq in FastqParser(input_file, verbose, trim_bad_bases=True, make_seq=make_seq):
        i += 1
        if len(seq) < minimum_length:
            num_too_small += 1
            continue

        if i % 100000 == 0:
            print name
            print i

        data = "\n".join(['>%s' % name, seq, ''])
        if not test_run:
            outfile.write(data)
        else:
            print data

        if test_run and i > 10000:
            break

    if not test_run:
        outfile.close()

    print '\nSuccess!!\n%d Sequences were read' % i
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
                    default=40,
                    help='minimum length of sequences to write [default: %default]')
                    ]

    parser, opts, args = parse_command_line_parameters(**script_info)
    run(opts.input_file, opts.output_file, opts.minimum_length, opts.test_run, opts.verbose)