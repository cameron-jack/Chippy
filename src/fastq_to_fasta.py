import cl_options as clo
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

def config_options(parser):
    """sets up the command line options"""
    
    parser.add_option('-i', '--input_file',
       help="name of the input file")
    
    parser.add_option('-o', '--output_file',
       help="name of the output results file")
    
    parser.add_option('-l', '--minimum_length', type='int', default=40,
        help = 'minimum length of sequences to write [default: %default]')
    
    parser.add_option('-t', '--test_run', action='store_true', default=False,
        help="test run, [default: %default]")

if __name__ == "__main__":
    opts, args = clo.parse_command_line_parameters(config_options,
                    required_options=['input_file', 'output_file'])
    run(opts.input_file, opts.output_file, opts.minimum_length, opts.test_run,
        opts.verbose)
