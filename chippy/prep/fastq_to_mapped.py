#!/usr/bin/env python
"""implements the full workflow of processing sequences
- converts to fasta in prep for using blat
- uses blat to find adapter sequences
- produces trimmed sequences without adpaters
- maps the clean sequences using bowtie
"""
import sys
sys.path.extend(['..', '../src'])

import numpy

from cogent.util.misc import parse_command_line_parameters
from optparse import make_option

from chippy.prep import reduce, pristine_seqs, command_line, fastq_to_fasta
from chippy.prep.mapped_files import MappedFiles
from chippy.util.run_record import RunRecord

__author__ = "Gavin Huttley, Cameron Jack"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley, Cameron Jack"
__credits__ = ["Gavin Huttley, Cameron Jack"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"
__version__ = '0.1'

script_info = fastq_to_fasta.script_info

optional = script_info['optional_options']
optional.append(make_option('-p', '--pval_cutoff',
                type='float', default=0.001,
                help='Minimum p-value for mapping quality of reads\
                [default: %default]'))

required = script_info['required_options']
required.append(make_option('--blat_adapters',
                            help='path to the Illumina adapters'))
required.append(make_option('--adapter_clipper',
                            choices=['blat','fastx'],
                            help='Choose adapter clipping program'))
required.append(make_option('--index',
                            help='path to the aligner genome index'))
required.append(make_option('--aligner', type='choice',
                            choices=['bwa', 'bowtie'],
                            help='Choose the aligner program'))

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
    mapped_files = MappedFiles(opts.input_file, opts.save_dir)
    fastq_fn = mapped_files.fastq_fn
    fasta_fn = mapped_files.fasta_fn
    working_dn = mapped_files.working_dn
    psl_fn = mapped_files.psl_fn
    trimmed_fn = mapped_files.trimmed_fn

    contaminated_fn = mapped_files.contaminated_fn
    pristine_fn =  mapped_files.pristine_fn
    combined_fn = mapped_files.combined_fn

    combined_sai_fn = mapped_files.combined_sai_fn
    combined_sam_fn = mapped_files.combined_sam_fn

    run_record_fn = mapped_files.run_record_fn

    # run_records tracks metrics from each step that are printed & saved at
    # completion of the entire workflow. Expected format is: program_name
    run_record = RunRecord()

    if opts.adapter_clipper.lower() == 'blat':
        # convert fastq to fasta so we can use blat to map the adapters
        # this function call creates the save_dir path, if it doesn't already
        # exist. Rewrite fastq=true
        print 'Prepping for blat (fastq to fasta)'
        run_record = fastq_to_fasta.run(fastq_fn, fasta_fn,
            opts.minimum_length, True, run_record, opts.test_run)
    
        print 'Running blat'
        run_record = command_line.run_blat(opts.blat_adapters,
                fasta_fn, psl_fn, run_record, opts.test_run)
    
        print 'Writing seqs without adapters'
        run_record = pristine_seqs.main(psl_fn, trimmed_fn, run_record,
                     opts.test_run)

        # concatenate the pristine and contaminated fastq files
        print 'Concatenating contaminated and pristine fastq files'
        run_record = command_line.concatenate(pristine_fn, contaminated_fn,
                    combined_fn, run_record, opts.test_run)
    elif opts.adapter_clipper.lower() == 'fastx':
        # We aren't separating pristine and contaminated with fastx
        # so write straight to 'combined' file
        print 'Running fastx adapter clipping'
        run_record = command_line.run_fastx_clipper(opts.blat_adapters,
                    fastq_fn, trimmed_fn, run_record, opts.test_run)

        print 'Running fastq_quality_trimmer'
        run_record = command_line.run_fastq_qual_trim(trimmed_fn, combined_fn,
                     run_record, opts.test_run)
    else:
        raise RuntimeError('Unknown adapter clipper choice %s' \
                % opts.adapter_clipper)

    if opts.aligner.lower() == 'bowtie':
        print 'Mapping seqs with bowtie'
        run_record = command_line.run_bowtie(opts.index,
            combined_fn, combined_sam_fn, run_record, opts.test_run)
    elif opts.aligner.lower() == 'bwa':
        print 'Aligning seqs with bwa aln'
        run_record = command_line.run_bwa_aln(opts.index,
            combined_fn, combined_sai_fn, run_record, opts.test_run)
        print 'Finding chromosomal coordinates with bwa samse'
        run_record = command_line.run_bwa_samse(opts.index, 
            combined_sai_fn, combined_fn, combined_sam_fn,
            run_record, opts.test_run)
    else:
        raise RuntimeError('Unknown aligner choice %s' % opts.aligner)

    # minimal_read map the concatenated file
    run_record = reduce.run(infile_name=combined_sam_fn,
        outdir=opts.save_dir, chroms='Do All', pval_cutoff=opts.pval_cutoff,
        limit=numpy.inf, run_record=run_record, dry_run=opts.test_run)
    # display/write synopsis of run
    run_record.display()
    table = run_record.getMessageTable()
    table.writeToFile(run_record_fn, sep='\t')

if __name__ == "__main__":
    main()

