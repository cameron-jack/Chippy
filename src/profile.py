#!/usr/bin/env python

import sys
sys.append('..')

try:
    import line_profiler
except ImportError:
    print 'Requires line_profiler module'
    print '$ sudo easy_install -U line_profiler'
    exit()

#from plot_quality import plot_quality, make_display
#from inline_stats import RunningStats
#from light_seq import LightSeq

from chippy.prep import pristine_seqs

psl_file = '../../../tremethick/blat_output/s_6_blat.psl'
fastq_file = '../../../tremethick/trimmed/s_6_sequence_trimmed.fastq'
output = 'test_profile'

profiler = line_profiler.LineProfiler(pristine_seqs.run, pristine_seqs.get_corrupt_seq_names, pristine_seqs.write_pristine)
profiler.run("pristine_seqs.run('%s', '%s', '%s', %s)"%(psl_file, fastq_file, output, 'False'))
profiler.print_stats()
