#!/usr/bin/env python
"""implements the full workflow of processing sequences
- converts to fasta in prep for using blat
- uses blat to find adapter sequences
- produces trimmed sequences without adpaters
- maps the clean sequences using bowtie
"""

import os
import time
import numpy

from cogent import LoadTable
from cogent.util.misc import parse_command_line_parameters
from optparse import make_option

from chippy.prep import reduce, pristine_seqs, command_line, fastq_to_fasta

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2011, Gavin Huttley"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"
__version__ = '0.1'

class RunRecord(object):
    """object for recording program messages"""
    def __init__(self):
        super(RunRecord, self).__init__()
        self.records = []
        
    def addMessage(self, program_name, error_type, message, value):
        """add a message about an execution"""
        self.records.append([program_name, error_type, message, value])
    
    def getMessageTable(self):
        """docstring for display"""
        header = ['program_name', 'message type', 'message', 'value']
        table = LoadTable(header=header, rows=self.records)
        return table
    
    def display(self):
        table = self.getMessageTable()
        print table


def make_fastq_output_filename(input_file):
    """makes the fasta output filename from the fastq input name"""
    return input_file.replace('_sequence.txt', '.fasta')

script_info = fastq_to_fasta.script_info
required = script_info['required_options']
required.append(make_option('--blat_adapters',
                            help='path to the Illumina adapters'))
required.append(make_option('--bowtie_index',
                            help='path to the bowtie genome index'))

# output_file is required for fastq_to_fasta, but not as part of full chain
# since we will automatically generate it
output_file = None
for option in required:
    if option.dest == 'output_file':
        output_file = option
        break

assert output_file is not None
required.remove(output_file)

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    
    # make the assorted filenames
    output_file = make_fastq_output_filename(opts.input_file)
    working_dir = '%s-working' % (opts.save_dir)
    fasta_file = os.path.join(working_dir, output_file)
    psl_out = fasta_file.replace('.fasta', '.psl')
    trimmed_fastq = fasta_file.replace('.fasta', '_trimmed.fastq')
    pristine_fastq = trimmed_fastq.replace('_trimmed',
                                    '_trimmed_pristine')
    pristine_map = pristine_fastq.replace('.fastq', '.map')
    contaminated_fastq = pristine_fastq.replace('pristine', 'contaminated')
    contaminated_map = contaminated_fastq.replace('.fastq', '.map')
    combined_map = fasta_file.replace('.fasta', '.map')
    run_record_file_name = os.path.join(opts.save_dir, 'run_record.txt')
    
    # run_records tracks metrics from each step that are printed & saved at
    # completion of the entire workflow. Expected format is: program_name
    run_record = RunRecord()
    # convert fastq to fasta so we can use blat to map the adapters
    # this function call creates the save_dir path, if it doesn't already
    # exist
    print 'Prepping for blat (fastq to fasta)'
    run_record = fastq_to_fasta.run(opts.input_file, working_dir,
        output_file, opts.minimum_length, True,
        run_record, opts.test_run)
    
    print 'Running blat'
    run_record = command_line.run_blat(opts.blat_adapters,
            fasta_file, psl_out, run_record, opts.test_run)
    
    print 'Writing seqs without adapters'
    run_record = pristine_seqs.main(psl_out, trimmed_fastq, run_record,
            opts.test_run)
    
    print 'Mapping contaminated seqs'
    run_record = command_line.run_bowtie(opts.bowtie_index, working_dir,
        contaminated_fastq, contaminated_map, run_record, opts.test_run)
    print 'Mapping pristine seqs'
    run_record = command_line.run_bowtie(opts.bowtie_index, working_dir,
        pristine_fastq, pristine_map, run_record, opts.test_run)
    
    # concatenate the pristine and contaminated files
    print 'Concatenating contaminated and pristine map files'
    run_record = command_line.concatenate(pristine_map, contaminated_map,
                combined_map, run_record, opts.test_run)
    
    # minimal_read map the concatenated file
    run_record = reduce.run(input_file=combined_map,
        outdir=opts.save_dir, chroms='Do All', limit=numpy.inf,
        run_record=run_record, dry_run=opts.test_run)
    # display/write synopsis of run
    run_record.display()
    table = run_record.getMessageTable()
    table.writeToFile(run_record_file_name, sep='\t')

if __name__ == "__main__":
    main()
