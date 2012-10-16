#!/usr/bin/env python

import sys, os
sys.path.append('..')

from shutil import rmtree

from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files

from chippy.prep.command_line import run_command

class FastqPipe(TestCase):
    def test_pipe(self):
        """Pipeline changes should produce same results as original pipeline"""

        # Check that we're not going to nuke something we shouldn't
        if os.path.exists('test_pipe_dir'):
            rmtree('test_pipe_dir')
        if os.path.exists('test_pipe_dir-working'):
            rmtree ('test_pipe_dir-working')

        # Run pipeline with known fastq file
        call_string = 'python2.7 ../chippy/prep/fastq_to_mapped.py -i data/known_sequence.fq.gz '\
                    '-S test_pipe_dir --blat_adapters $ADAPTERS '\
                    '--adapter_clipper fastx --index $BWA_62_mouse_INDEX '\
                    '--aligner bwa' # -t flag not fully implemented yet, raises a different bug

        returncode, stdout, stderr = run_command(call_string)

        if returncode != 0:
            raise RuntimeError('fastq_to_mapped exited with error: ' + str(stderr))
        # Generate md5 checksums of result files
        call_string = 'gunzip -f test_pipe_dir/*.gz'

        returncode, stdout, stderr = run_command(call_string)
        if returncode != 0:
            raise RuntimeError('gunzip exited with an error: ' + str(stderr))

        call_string = 'md5sum test_pipe_dir/*.txt > test_pipe_dir/results.md5'
        
        returncode, stdout, stderr = run_command(call_string)
        if returncode != 0:
            raise RuntimeError('md5sum exited with an error: ' + str(stderr))


        # Compare values. Note: this could fail with user altered parameters
        resulting_md5 = open ('test_pipe_dir/results.md5')
        expected_md5 = open ('data/expected_seq_results.md5')
        for line_e, line_r in zip(resulting_md5, expected_md5):
            self.assertEqual(line_e, line_r)

        resulting_md5.close()
        expected_md5.close()
        
        # now clean up dirs
        rmtree ('test_pipe_dir')
        rmtree ('test_pipe_dir-working')

if __name__ == "__main__":
    main()
