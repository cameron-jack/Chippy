#!/usr/bin/env python
"""implements the full workflow of processing sequences
- converts to fasta in prep for using blat
- uses blat to find adapter sequences
- produces trimmed sequences without adpaters
- maps the clean sequences using bowtie
"""
import sys, os, shutil, zipfile, bz2, gzip
sys.path.extend(['..', '../src'])

from cogent.util.misc import parse_command_line_parameters
from optparse import make_option

from chippy.prep import command_line
from chippy.prep.mapped_files import MappedFiles, MappedSangerFiles
from chippy.util.run_record import RunRecord

__author__ = "Gavin Huttley, Aaron, Chuah, Cameron Jack"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley, Cameron Jack, Aaron Chuah"
__credits__ = ["Gavin Huttley, Cameron Jack, Aaron Chuah"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"
__version__ = '0.1.3'

script_info = {}
descr = "Draft snp discovery process"
script_info['brief_description']= descr
script_info['script_description'] = descr
script_info['version'] = '1e-3.alpha'
script_info['script_usage']=[]


script_info['help_on_no_arguments'] = True
script_info['required_options'] = [
    make_option('-1','--input_file_1',
                help='The 1st input fastq sequence file'),
    make_option('-2','--input_file_2',
                help='The 2nd input fastq sequence file'),
    make_option('-S', '--save_dir', help='path to save all files'),
    make_option('--adapters', help='path to the Illumina adapters'),
    make_option('--index', help='path to the bwa aligner index')
]

script_info['optional_options'] = [\
    make_option('-n', '--num_threads', type='int',
                default = 6,
                help='Number of threads to use [default: %default]'),
    make_option('-m', '--mem_usage', type='int',
                default = 7000000000,
                help='memory usage for samtools sort [default: %default]'),
    make_option('-t', '--test_run', action='store_true',
                dest= 'test_run', default = False,
                help='Dry run without writing any data'
                +' [default: %default]'),
    make_option('-s', '--sample_name', type='string', default = '',
                help='specify sample name in annotated bam header'
                +' [default: %default]'),
    make_option('-w', '--work_dir', type='string', default = '',
                help='specify temporary working directory'
                +' [default: %default]'),
    make_option('-b', '--begin', type='int', default = 1,
                help='begin at stage # [default: %default]'),
    make_option('-e', '--end', type='int', default = 7,
                help='end at stage # [default: %default]'),
    make_option('-D', '--delete', action='store_true', default = False,
                help='Deletes the working dir at the end of the run'
                +' [default: %default]'),
    make_option('-r', '--reduce', action='store_true', default = False,
                help='Finish with ChIP-Seq reduction step to old format'
                +' [default: %default]'),
    make_option('-I', '--Illumina_version', default=1.7,
                help='Illumina pipeline version number'
                +' [default: %default]'),
    make_option('-p', '--pval_cutoff', type='float', default=0.05,
        help='Minimum p-value for mapping quality of reads '\
             '[default: %default]')
]

def main():
    option_parser, opts, args =\
        parse_command_line_parameters(**script_info)

    rr = RunRecord() # create the run audit tracker

    ### We support different two very different pipelines depending on Illumina version

    # Make the assorted filenames required

    ## v1.3 - v1.7 pipeline:
    # Zipped or unzipped, remove 5' adapters with fastx_clipper on each end
    # Remove poor quality bases and reads with fastx_quality_trimmer on each end
    # Balance unmatched read pairs
    # Run bwa aln on each end
    # Run bwa sampe and save as either BAM or reduced

    if 1.3 >= opts.Illumina_version <= 1.7:
        mapped_files_1 = MappedFiles(opts.input_file_1, opts.save_dir, opts.work_dir)
        filenames_1 = dict(fastq=mapped_files_1.adapterless_trimmed_fn,
                       pristine=mapped_files_1.pristine_fn, sai=mapped_files_1.sai_fn,
                       bam=mapped_files_1.bam_fn, bed=mapped_files_1.bed_fn,
                        work_dir=mapped_files_1.working_dn)
        mapped_files_2 = MappedFiles(opts.input_file_2, opts.save_dir, opts.work_dir)
        filenames_2 = dict(fastq=mapped_files_2.adapterless_trimmed_fn,
                       pristine=mapped_files_2.pristine_fn, sai=mapped_files_2.sai_fn,
                       bam=mapped_files_2.bam_fn)

    ## v1.8+
    # Check if gzipped, if yes then unzip - raw
    # Run Perl DynamicTrim.pl on each end (outputs as filename.ext.trimmed) - trimmed
    # Run Sickle on both ends (we name the out files) - pristine
    # Run bwa aln on each end - sai
    # Run bwa sampe and save as either BAM or reduced

    elif opts.Illumina_version >= 1.8: # Illumina pipeline 1.8 or higher
        mapped_files_1 = MappedSangerFiles(opts.input_file_1, opts.save_dir, opts.work_dir)
        filenames_1 = dict(unzipped=mapped_files_1.unzipped_fq_fn,
                trimmed=mapped_files_1.trimmed_fn,
                pristine=mapped_files_1.pristine_fn, sai=mapped_files_1.sai_fn)

        mapped_files_2 = MappedSangerFiles(opts.input_file_2, opts.save_dir, opts.work_dir)
        filenames_2 = dict(unzipped=mapped_files_2.unzipped_fq_fn,
                trimmed=mapped_files_2.trimmed_fn,
                pristine=mapped_files_2.pristine_fn, sai=mapped_files_2.sai_fn)

        # unspecific to an 'end' - they're still derived from mapped_files_1
        unified_fns = dict(unsorted=mapped_files_1.unsorted_bam_fn,
                sorted=mapped_files_1.sorted_bam_fn,
                filtered=mapped_files_1.filtered_bam_fn,
                bed=mapped_files_1.bed_fn,
                work_dir=mapped_files_1.working_dn,
                run_record=mapped_files_1.run_record_fn)

    else: # unsupported Illumina pipeline version
        rr.addError('fastq_to_mapped_bam', 'Invalid Illumina pipeline given. v1.3+ supported',
            opts.Illumina_version)
        rr.display()
        sys.exit(0)

    # Make sure we don't overwrite any existing output BAM file
    if not opts.reduce and os.path.exists(filenames_1['bam']):
            rr.addInfo('fastq_to_mapped_bam', 'already exists, exiting', filenames_1['bam'])
            rr.display()
            sys.exit(0)

    # Quality control pipeline - depending on input data version
    if 1.3 >= opts.Illumina_version <= 1.7:

        # do input file 1
        if opts.begin <= 1 and opts.end >= 1:
            rr = command_line.run_fastx_clip_and_trim(opts.adapters, opts.input_file_1,
                        filenames_1['fastq'], rr, opts.num_threads, opts.test_run)

        # do input file 2
        if opts.begin <= 2 and opts.end >= 2:
            rr = command_line.run_fastx_clip_and_trim(opts.adapters, opts.input_file_2,
                        filenames_2['fastq'], rr, opts.num_threads, opts.test_run)

        # produce the paired pristine seq files with matched reads
        if opts.begin <= 3 and opts.end >= 3:
            # rr = command_line.run_pristine_paired(
            rr = command_line.run_pristine_seesaw(opts.input_file_1,
                        filenames_1['fastq'], filenames_2['fastq'],
                        filenames_1['pristine'], filenames_2['pristine'],
                        rr, opts.num_threads, opts.test_run)

    elif opts.Illumina_version >= 1.8: # Illumina pipeline 1.8 or higher
        # Check if gzipped, bzipped or zipped, and if yes, then unzip to work directory
                # otherwise copy so DynamicTrim can work in-place
        if opts.begin <= 1 and opts.end >= 1:
            # input file 1
            if mapped_files_1.is_gzip:
                fastq_out = open(filenames_1['unzipped'], 'w')
                for line in gzip.GzipFile(opts.input_file_1, 'rb'):
                    fastq_out.write(line)
                fastq_out.close()
                rr.addInfo('fastq_to_mapped_bam', 'file decompressed with gz to work dir', opts.input_file_1)
            elif mapped_files_1.is_bzip2:
                fastq_out = open(filenames_1['unzipped'], 'w')
                for line in bz2.BZ2File(opts.input_file_1, 'rb'):
                    fastq_out.write(line)
                fastq_out.close()
                rr.addInfo('fastq_to_mapped_bam', 'file decompressed with bz2 to work dir', opts.input_file_1)
            # Zip support seems broken right now, so removing.
            #elif mapped_files_1.is_zip:
            #    zip1 = zipfile.ZipFile(opts.input_file_1, 'r')
            #    zip1.extractall(opts.work_dir)
            #    rr.addInfo('fastq_to_mapped_bam', 'file unzipped to work dir', opts.input_file_1)
            elif mapped_files_1.is_fastq:
                shutil.copy(opts.input_file_1, filenames_1['unzipped'])
                rr.addInfo('fastq_to_mapped_bam', 'file copied to work dir', opts.input_file_1)
            else:
                rr.addError('fastq_to_mapped_bam', 'unrecognised filename extension', opts.input_file_1)

            # input file 2
            if mapped_files_2.is_gzip:
                fastq_out = open(filenames_2['unzipped'], 'w')
                for line in gzip.GzipFile(opts.input_file_2, 'rb'):
                    fastq_out.write(line)
                fastq_out.close()
                rr.addInfo('fastq_to_mapped_bam', 'file decompressed with gz to work dir', opts.input_file_2)
            elif mapped_files_2.is_bzip2:
                fastq_out = open(filenames_2['unzipped'], 'w')
                for line in bz2.BZ2File(opts.input_file_2, 'rb'):
                    fastq_out.write(line)
                fastq_out.close()
                rr.addInfo('fastq_to_mapped_bam', 'file decompressed with bz2 to work dir', opts.input_file_2)
            # Zip support seems broken right now, so removing
            #elif mapped_files_2.is_zip:
            #    zip2 = zipfile.ZipFile(opts.input_file_2, 'r')
            #    zip2.extractall(opts.work_dir)
            #    rr.addInfo('fastq_to_mapped_bam', 'file unzipped to work dir', opts.input_file_2)
            elif mapped_files_2.is_fastq:
                shutil.copy(opts.input_file_2, filenames_2['unzipped'])
                rr.addInfo('fastq_to_mapped_bam', 'file copied to work dir', opts.input_file_2)
            else:
                rr.addError('fastq_to_mapped_bam', 'unrecognised filename extension', opts.input_file_2)

        if opts.begin <= 2 and opts.end >= 2:
            # Don't worry about 5' adapter clipping for Hi-Seq, though we might want to think about 3' adapters in future

            # Remove low quality bases with DynamicTrim (SolexaQA)
            rr = command_line.run_dynamic_trim(filenames_1['unzipped'], unified_fns['work_dir'],
                    opts.pval_cutoff, rr, opts.test_run)
            rr = command_line.run_dynamic_trim(filenames_2['unzipped'], unified_fns['work_dir'],
                    opts.pval_cutoff, rr, opts.test_run)

        if opts.begin <= 3 and opts.end >= 3:
            # Balance inputs files with Sickle (UC Davis), discard single-ends residual
            min_length = 19;
            rr = command_line.run_sickle_pe(filenames_1['trimmed'], filenames_2['trimmed'],
                        filenames_1['pristine'], filenames_2['pristine'],
                        min_length, rr, opts.test_run)

    else:
        rr.addError('fastq_to_mapped_bam', 'Invalid Illumina pipeline given. v1.3+ supported',
                    opts.Illumina_version)
        rr.display()
        sys.exit(0)

    # Now align each individually
    if opts.begin <= 4 and opts.end >= 4:
        rr = command_line.run_bwa_aln(opts.index, filenames_1['pristine'],
                filenames_1['sai'], opts.Illumina_version,
                rr, opts.num_threads, opts.test_run)

    if opts.begin <= 5 and opts.end >= 5:
        rr = command_line.run_bwa_aln(opts.index, filenames_2['pristine'],
                filenames_2['sai'], opts.Illumina_version,
                rr, opts.num_threads, opts.test_run)

    # Note - from here the pipeline is broken for Illumina version < 1.8
    if opts.begin <= 6 and opts.end >= 6:
        rr = command_line.bwa_sampe_to_bam_with_header(opts.index,
                filenames_1['sai'], filenames_2['sai'],
                filenames_1['pristine'], filenames_2['pristine'],
                unified_fns['unsorted'],
                rr, opts.input_file_1, opts.test_run)

    # Sort BAM
    if opts.begin <= 7 and opts.end >= 7:
        rr = command_line.bam_to_sorted_bam(unified_fns['unsorted'],
                unified_fns['sorted'], opts.mem_usage, rr, opts.test_run)

    ## index bam or filter and output to BED for ChIP-Seq analysis
    if opts.begin <= 8 and opts.end >= 8:
        if opts.reduce:
            # Reduce.py is now legacy for single-end reads

            # filter out bad mappings and PCR duplicates
            rr = command_line.convert_sorted_bam_to_filtered_bam(
                    unified_fns['sorted'], unified_fns['filtered'], rr,
                    opts.test_run)

            # convert BAM to BED
            rr = command_line.convert_bam_to_bed(unified_fns['filtered'],
                    unified_fns['bed'], rr, opts.test_run)

            # Single command version of the above
            #rr = command_line.convert_sorted_bam_to_filtered_bed(filenames_1['bam'],
            #        filenames_1['bed'], rr, opts.test_run)

        else:
            rr = command_line.index_bam(filenames_1['bam'], rr, opts.test_run)


    if opts.delete:
        shutil.rmtree(filenames_1['work_dir'])

    ## output audit
    rr.display()
    table = rr.getMessageTable()
    # print "\n\nRun record location: %s\n\n" % unified_fns['run_record']
    table.writeToFile(unified_fns['run_record'], sep='\t')

if __name__ == "__main__":
    main()