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

class MappedFiles:

    # Self to pointer to single instance
    __single = None

    def make_fasta_from_fastq(self, fastq_filename):
        """makes the fasta output filename from the fastq input name"""
        fastq_filename = os.path.basename(fastq_filename)
        pattern = '(_sequence\.txt|\.fastq|\.fq)$'
        if not re.search(pattern, fastq_filename):
            raise RuntimeError(
            "Input file name [%s] doesn't match expected convention" % fastq_filename)

        fasta_filename = re.sub(pattern, '.fasta', fastq_filename)
        return fasta_filename

    def __init__(self, fastq_fname, save_dname):
        """These files need the following filenames:
        fastq_to_mapped: input.fastq
        """
        ### These 3 code lines enforce the Singleton pattern
        if MappedFiles.__single:
            raise MappedFiles.__single

        MappedFiles.__single = self

        # File names needed globally:
        self.fastq_fn = fastq_fname
        fasta_fn = self.make_fasta_from_fastq(fastq_fname)

        # directory names, make sure they exist
        create_path(save_dname)
        self.save_dn = save_dname
        self.working_dn = '%s-working' % save_dname
        create_path(self.working_dn)

        # All file names from here include the working path:
        self.fasta_fn = os.path.join(self.working_dn, fasta_fn)
        self.psl_fn = self.fasta_fn.replace('.fasta', '.psl')
        self.trimmed_fn = self.fasta_fn.replace('.fasta',
                                           '_trimmed.fastq')

        self.contaminated_fn = self.fasta_fn.replace('.fasta',
                                                '_contaminated.fastq')
        self.pristine_fn = self.fasta_fn.replace('.fasta',
                                            '_pristine.fastq')
        self.combined_fn = self.fasta_fn.replace('.fasta', '_combined.fastq')
        self.combined_sai_fn = self.combined_fn.replace('.fastq', '.sai')
        self.combined_sam_fn = self.combined_fn.replace('.fastq', '.sam')

        self.run_record_fn = os.path.join(save_dname, 'run_record.txt')

def mapped_file_handle(fastq_fname, save_dname):
    """mapped_file_handle provides the safe way of creating an instance of MappedFiles
    or returning THE existing instance of MappedFiles"""
    try:
        singleton = MappedFiles(fastq_fname, save_dname)
    except MappedFiles, m:
        singleton = m
    return singleton

def test_instantiation():
    """used for testing distance vs local instantiation and for returning
    an instance after initial creation"""
    mapped_files = mapped_file_handle(TEST_FILE_NAME, TEST_DIR_NAME)
    return mapped_files

def internal_instantiation():
    """used for testing distance vs local instantiation and for returniTEST_FILE_NAMEng
    an instance after initial creation"""
    mapped_files = mapped_file_handle(TEST_FILE_NAME, TEST_DIR_NAME)
    if mapped_files.fastq_fn == TEST_FILE_NAME:   
        raise RuntimeError("mapped_files not instantiated. \
        Use mapped_file_handle(filename, dirname)")
    return mapped_files

def main():
    """A self test of both being a singleton and also of the usage and outputs """
    
    mapped_files = mapped_file_handle("inputFile.fastq", "/blah")

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

    mapped_files = mapped_file_handle("inputFile2.fastq", "/blah2")

    print mapped_files.fastq_fn

if __name__ == "__main__":
    main()