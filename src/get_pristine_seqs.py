from parse_psl import MinimalPslParser
from parse_fastq import MinimalFastqParser
from light_seq import LightSeq

def get_corrupt_seq_names(psl_name, test_run):
    psl_parser = MinimalPslParser(psl_name)
    psl_parser.next()
    psl_parser.next()

    num = 0
    contaminated_info = []
    for record in psl_parser:
        # get the contaminated sequence name and the index where contamination
        # starts
        info = record[9], record[11]
        num += 1
        contaminated_info.append(info)
        if test_run and num >= 1000:
            break

    return contaminated_info

to_fasta = lambda name, seq: '\n'.join(['>%s' % name, seq])

def write_pristine(fastq_name, outfile_root, not_pristine, test_run):
    num = 0
    if not test_run:
        outfile_pristine = open(outfile_root + '_pristine.fastq', 'w')
        outfile_contaminated = open(outfile_root + '_contaminated.fastq', 'w')

    not_pristine_names = [a[0] for a in not_pristine]

    for name, seq, qual in MinimalFastqParser(open(fastq_name)):
        num += 1
        seq_object = LightSeq(seq, name, qual)


        if test_run:
            print fastq_formatted
            if num > 100:
                break

        # if sequence is contaminated with adapater, trim the sequence
        if name in not_pristine_names:
            index = not_pristine_names.index(name)
            start = not_pristine[index][1]

            # since we trim everything after the adapter sequence start,
            # there is little value in keeping a sequence which is going to be
            # less than 35 bp long.
            if start < 35:
                continue

            seq_object = seq_object[:start]
            fastq_formatted = seq_object.toFastq()
            outfile_contaminated.write(fastq_formatted + '\n')

        else:
            fastq_formatted = seq_object.toFastq()
            outfile_pristine.write(fastq_formatted + '\n')

    if not test_run:
        outfile_contaminated.close()
        outfile_pristine.close()

def run(input_psl_file, input_file, outfile_root, test_run):
    """identifies reads not to be written, then writes everything else"""
    not_pristine = get_corrupt_seq_names(input_psl_file, test_run)
    write_pristine(input_file, outfile_root, not_pristine, test_run)
    print '\n Done!'

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
        """python get_pristine_seqs.py -p <somefile.psl> -i <seqs> -o <outfile_root> -t"""))

    script_info['help_on_no_arguments'] = True
    script_info['required_options'] = [
        make_option('-p','--input_psl_file',
                    help='The input psl file from blat'),
        make_option('-i','--input_file',
                    help='The input fasta sequence file'),
        make_option('-o','--outfile_root',
                    help='The pristine and contaminated files will be written '\
                    'to two different files, with this as this + _pristine or '\
                    '_contaminated as the name of the file'),
        ]

    script_info['optional_options'] = [\
        make_option('-t','--test_run', action='store_true',
                    dest='test_run', default = False,
                    help='Dry run without writing any data'
                    +'[default: %default]'),
                    ]

    parser, opts, args = parse_command_line_parameters(**script_info)

    run(opts.input_psl_file, opts.input_file, opts.outfile_root, opts.test_run)
