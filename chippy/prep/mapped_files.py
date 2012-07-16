#!/usr/bin/env python
"""MappedFiles: returns filename strings based on standard object names.
Implemented using the Singleton pattern described here:
http://www.python.org/workshops/1997-10/proceedings/savikko.html
"""

import os, re
from chippy.util.util import create_path

__author__ = "Cameron Jack"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley, Cameron Jack"
__credits__ = ["Cameron Jack"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"
__version__ = '0.1'

TEST_FILE_NAME = "testFile1.fastq"
TEST_NAME = "testFile1"
TEST_DIR_NAME = "/blah1"

# .fq.gz, .fastq.gz, .fq, .fastq

def make_fasta_from_fastq(fastq_filename):
    """makes the fasta output filename from the fastq input name"""
    fastq_filename = os.path.basename(fastq_filename)
    pattern = '(_sequence\.txt|\.fastq|\.fq)(\.gz)*$'
    
    fasta_filename = re.sub(pattern, '.fasta', fastq_filename)
    if fastq_filename.endswith('.gz'):
        fasta_filename += '.gz'
    
    return fasta_filename

def make_unified_fn(fastq_filename):
    """returns a filename with removal of paired end component"""
    uni = re.sub(r'_(1|2)(\.|_)', r'\2', fastq_filename)
    return uni

def robust_replace_suffix(fname, search, replace):
    if fname.endswith('.gz') and not search.endswith('.gz'):
        search = search + '.gz'
        if not replace.endswith('.gz'):
            replace = replace + '.gz'
    return fname.replace(search, replace)

class MappedFiles:


    def __init__(self, fastq_fname, save_dname, working_dname):
        """These files need the following filenames:
        fastq_to_mapped: input.fastq
        """
        if save_dname.endswith('/'):
            save_dname = save_dname[:-1]

        if working_dname.endswith('/'):
            working_dname = working_dname[:-1]

        if working_dname=='':
            working_dname = save_dname+'-working'

        # File names needed globally:
        self.fastq_fn = fastq_fname
        self.fasta_fn = os.path.join(working_dname, make_fasta_from_fastq(fastq_fname))
        self.psl_fn = self.fasta_fn.replace('.fasta', '.psl')
        self.contaminated_fn = self.fasta_fn.replace('.fasta', '_contaminated.fq')
        self.trimmed_fn = robust_replace_suffix(self.fasta_fn, '.fasta', '_trimmed.fq')
        fastq_path = fastq_fname.split('/')
        self.fastq_fn_only = fastq_path[-1]
        # directory names, make sure they exist
        create_path(save_dname)
        self.save_dn = save_dname
        self.working_dn = working_dname
        create_path(self.working_dn)

        # All file names from here include the working path:
        suffix = ['.fastq', '.fq']['.fq' in self.fastq_fn_only]
        
        self.adapterless_trimmed_fn = self.working_dn+'/'+robust_replace_suffix(self.fastq_fn_only, suffix,
            '_adapterless_trimmed.fq')
        ###
        self.pristine_fn = robust_replace_suffix(self.adapterless_trimmed_fn, '_adapterless_trimmed.fq', '_pristine.fq')
        #suffix = ['.fq', '.fq.gz'][self.pristine_fn.endswith('.gz')]
        self.sai_fn = robust_replace_suffix(self.pristine_fn, '_pristine.fq', '.sai')
        self.bam_fn = make_unified_fn(self.sai_fn.replace('.sai', '.bam')).replace(self.working_dn, self.save_dn)
        self.sam_fn = make_unified_fn(self.sai_fn.replace('.sai', '.sam')).replace(self.working_dn, self.save_dn)

        self.run_record_fn = self.bam_fn.replace('.bam','.run_record.txt')
        
        self.combined_fn = self.fasta_fn.replace('.fasta', '_combined.fq')
        self.combined_sai_fn = self.combined_fn.replace('.fq', '.sai')
        self.combined_sam_fn = self.combined_fn.replace('.fq', '.sam')

class MappedSangerFiles:

    def __init__(self, input_fastq_fname, save_dname, working_dname):
        """Create file and dir name strings for paired-end fastq files with
            Sanger scores"""

        if save_dname.endswith('/'):
            save_dname = save_dname[:-1]

        if working_dname.endswith('/'):
            working_dname = working_dname[:-1]

        if working_dname=='':
            working_dname = save_dname+'-working'

        # directory names, make sure they exist
        create_path(save_dname)
        self.save_dn = save_dname
        self.working_dn = working_dname
        create_path(self.working_dn)

        # The work flow is as follows:
        # Check if gzipped, if yes then unzip
        # Run Perl DynamicTrim.pl on each end (outputs as filename.ext.trimmed)
        # Run Sickle on both ends (we name the out files plus the singles file)
        # Run bwa aln on each end
        # Run bwa sampe and save bam
        # Sort bam - unsorted to sorted
        # Filter bam - sorted to filtered
        # Dedup bam - filtered to deduplicated
        # Convert BAM to BED (for ChipPy)

        # keep the full path in case needed throughout code
        self.input_fastq_fullpath = input_fastq_fname

        fastq_path = input_fastq_fname.split('/')
        # mostly we'll use the file name without the path
        fastq_fn_only = fastq_path[-1]

        # Set input file type
        self.is_gzip = False
        self.is_bzip2 = False
        self.is_zip = False # Code path not currently working
        self.is_fastq = False

        if fastq_fn_only.endswith('.gz'):
            fastq_fn_only_no_zip = fastq_fn_only.rstrip('.gz')
            self.is_gzip = True
        elif fastq_fn_only.endswith('.zip'):
            fastq_fn_only_no_zip = fastq_fn_only.rstrip('.zip')
            self.is_zip = True
        elif fastq_fn_only.endswith('.bz2'):
            fastq_fn_only_no_zip = fastq_fn_only.rstrip('.bz2')
            self.is_bzip2 = True
        elif fastq_fn_only.endswith('.fastq'):
            fastq_fn_only_no_zip = fastq_fn_only
            self.is_fastq = True
        elif fastq_fn_only.endswith('.fq'):
            fastq_fn_only_no_zip = fastq_fn_only
            self.is_fastq = True


        ## Working_dir names below

        self.unzipped_fq_fn = self.working_dn + '/' + fastq_fn_only_no_zip

        # trimmed is the result of DynamicTrim
        self.trimmed_fn = self.unzipped_fq_fn + '.trimmed'
        # balanced is the result of Sickle
        self.pristine_fn = self.unzipped_fq_fn + '.pristine'
        # singles is the remainder output of Sickle
        self.singles_fn = self.unzipped_fq_fn + '.singles'
        if self.unzipped_fq_fn.endswith('.fq'):
            # .sai is the bwa index file
            self.sai_fn = self.unzipped_fq_fn.replace('.fq', '.sai')
            # .unsorted.bam from bwa sampe
            self.unsorted_bam_fn = make_unified_fn(
                    self.unzipped_fq_fn.replace('.fq', '.unsorted.bam'))
            # .filtered.bam for ChipPy - keep just good quality paired-end mappings
            self.filtered_bam_fn = make_unified_fn(
                    self.unzipped_fq_fn.replace('.fq', '.filtered.bam'))
            # .filtered.bam for ChipPy - as filtered.bam but with PCR deduplication of reads performed.
            self.filtered_dedup_bam_fn = make_unified_fn(
                    self.unzipped_fq_fn.replace('.fq', '.filtered_dedup.bam'))
        elif self.unzipped_fq_fn.endswith('.fastq'):
            # .sai is the bwa index file
            self.sai_fn = self.unzipped_fq_fn.replace('.fastq', '.sai')
            # .unsorted.bam from bwa sampe
            self.unsorted_bam_fn = make_unified_fn(
                    self.unzipped_fq_fn.replace('.fastq', '.unsorted.bam'))
            # .filtered.bam - keep just good quality paired-end mappings
            self.filtered_bam_fn = make_unified_fn(
                    self.unzipped_fq_fn.replace('.fastq', '.filtered.bam'))
            # .filtered_dedup.bam - as filtered.bam but with PCR deduplication of reads performed.
            self.filtered_dedup_bam_fn = make_unified_fn(
                    self.unzipped_fq_fn.replace('.fastq', '.filtered_dedup.bam'))
        else:
            print 'In mapped_sanger_files, unrecognised filename extension '\
                    'for %s' % self.unzipped_fq_fn
            exit(0)

        ## Save_dir names below

        if self.unzipped_fq_fn.endswith('.fq'):
            # .sorted.bam is to be kept for all projects
            self.sorted_bam_fn = make_unified_fn(self.save_dn + '/' +\
                    (fastq_fn_only_no_zip.replace('.fq', '.sorted.bam')))
        elif self.unzipped_fq_fn.endswith('.fastq'):
            # .sorted.bam is to be kept for all projects
            self.sorted_bam_fn = make_unified_fn(self.save_dn + '/' +\
                    (fastq_fn_only_no_zip.replace('.fastq', '.sorted.bam')))

        # .sam is the final bwa output in tab delimited format - DEPRECATED
        self.sam_fn = self.sorted_bam_fn.replace('.sorted.bam','.sam')
        # .bed is a format needed for ChIP-Seq mappers and is now the standard format for ChipPy
        self.bed_fn = self.sorted_bam_fn.replace('.sorted.bam','.bed')
        # .gz.bed should be standard output for ChIP-Seq as it saves a lot of space and IO time
        # NOT CURRENTLY USED
        self.bed_gz_fn = self.sorted_bam_fn.replace('.sorted.bam', '.gz.bed')

        # output run_record to save directory
        self.run_record_fn = self.sorted_bam_fn.replace('.sorted.bam', '.run_record.txt')


def mapped_file_handle(fastq_fname, save_dname, work_dname):
    """mapped_file_handle provides the safe way of creating an instance of
    MappedFiles or returning THE existing instance of MappedFiles"""
    try:
        singleton = MappedFiles(fastq_fname, save_dname, work_dname)
    except MappedFiles, m:
        singleton = m
    return singleton

def test_instantiation():
    """used for testing distance vs local instantiation and for returning
    an instance after initial creation"""
    mapped_files = mapped_file_handle(TEST_FILE_NAME, TEST_DIR_NAME, '')
    return mapped_files

def internal_instantiation():
    """used for testing distance vs local instantiation and for returning
    an instance after initial creation"""
    mapped_files = mapped_file_handle(TEST_FILE_NAME, TEST_DIR_NAME, '')
    if mapped_files.fastq_fn == TEST_FILE_NAME:   
        raise RuntimeError("mapped_files not instantiated. \
        Use mapped_file_handle(filename, dirname)")
    return mapped_files

def main():
    """A self test of both being a singleton and also of the usage and outputs """
    
    mapped_files = mapped_file_handle("inputFile.fq", "/blah", "")

    print mapped_files.fastq_fn
    print mapped_files.fasta_fn
    print mapped_files.working_dn
    print mapped_files.psl_fn
    print mapped_files.contaminated_fn
    print mapped_files.pristine_fn
    print mapped_files.combined_fn
    print mapped_files.combined_sai_fn
    print mapped_files.combined_sam_fn
    print mapped_files.run_record_fn

    mapped_files = mapped_file_handle("inputFile2.fq", "/blah2", "")

    print mapped_files.fastq_fn

if __name__ == "__main__":
    # main()
    mapped_files = MappedFiles('/home/depressed/data/depression/HKC10038_HUMompR/CLEAN_FQ/PPRG_1045/test_1.bum.fq.gz',
        '/home/depressed/delme', "")
    print mapped_files.fastq_fn
    print mapped_files.adapterless_fn
    print mapped_files.adapterless_trimmed_fn
