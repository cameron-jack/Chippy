from cogent.parse.fasta import MinimalFastaParser

from parse_psl import MinimalPslParser

def get_corrupt_seq_names(psl_name, test_run):
    psl_parser = MinimalPslParser(psl_name)
    psl_parser.next()
    psl_parser.next()
    
    num = 0
    to_trim = set()
    for record in psl_parser:
        name = record[9]
        num += 1
        to_trim.update([name])
        if test_run and num >= 1000:
            break
    
    return to_trim

to_fasta = lambda name, seq: '\n'.join(['>%s' % name, seq])

def write_pristine(fasta_name, output_file, not_pristine, test_run):
    num = 0
    if not test_run:
        outfile = open(output_file, 'w')
    
    for name, seq in MinimalFastaParser(fasta_name):
        num += 1
        if name in not_pristine:
            continue
        
        fasta_formatted = to_fasta(name, seq)
        if test_run:
            print fasta_formatted
            if num > 100:
                break
        else:
            outfile.write(fasta_formatted + '\n')
    
    if not test_run:
        outfile.close()

def run(input_psl_file, input_file, output_file, test_run):
    """identifies reads not to be written, then writes everything else"""
    not_pristine = get_corrupt_seq_names(input_psl_file, test_run)
    write_pristine(input_file, output_file, not_pristine, test_run)
    print '\n Done!'

def _config_cl_parser(parser):
    """configures the command line options"""
    parser.add_option('-p', '--input_psl_file',
       help="name of the input psl file")
    
    parser.add_option('-i', '--input_file',
       help="name of the input fasta file")
    
    parser.add_option('-o', '--output_file',
       help="name of the output fasta file")
    
    parser.add_option('-t', '--test_run', action='store_true', default=False,
        help="test run, [default: %default]")
    

if __name__ == "__main__":
    from cl_options import parse_command_line_parameters
    
    opts, args = parse_command_line_parameters(_config_cl_parser,
        required_options=['input_psl_file', 'input_file', 'output_file'])
    
    run(opts.input_psl_file, opts.input_file, opts.output_file, opts.test_run)
