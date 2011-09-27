#!/usr/bin/env python

import sys, os
sys.path.append('..')

from cogent.util.unit_test import TestCase, main

from chippy.prep.command_line import run_command

class FastqPipe(TestCase):
    def test_pipe(self):
        """Pipeline changes should produce same results as original pipeline"""

        # Check that we're not going to nuke something we shouldn't
        if os.path.exists('test_pipe_dir'):
            raise RuntimeError('Cannot run test, "test_pipe_dir" already exists')

        # Run pipeline with known fastq file
        call_string = '../chippy/prep/fastq_to_mapped.py -i data/known_sequence.fq '\
                    '-S test_pipe_dir --blat_adapters $ADAPTERS '\
                    '--adapter_clipper fastx --index $BWA_62_mouse_INDEX '\
                    '--aligner bwa'

        returncode, stdout, stderr = run_command(call_string)

        # Generate md5 checksums of result files
        call_string = 'md5sum test_pipe_dir *.gz > test_pipe_dir/results.md5'
        returncode, stdout, stderr = run_command(call_string)

        # Compare values. Note: this could fail with user altered parameters
        resulting_md5 = open ('test_pipe_dir/results.md5')
        expected_md5 = open ('data/expected_seq_results.md5')
        for line_e, line_r in zip(resulting_md5, expected_md5):
            self.assertEqual(line_e, line_r)

        resulting_md5.close()
        expected_md5.close()
        
        # now clean up dirs
        call_string = 'rm -rf test_pipe_dir*'
        returncode, stdout, stderr = run_command(call_string)

if __name__ == "__main__":
    main()
