import sys
try:
    import line_profiler
except ImportError:
    print 'Requires line_profiler module'
    print '$ sudo easy_install -U line_profiler'
    exit()

#from plot_quality import plot_quality, make_display
#from inline_stats import RunningStats
#from light_seq import LightSeq

import get_pristine_seqs

psl_file = '../../../tremethick/blat_output/s_6_blat.psl'
fastq_file = '../../../tremethick/trimmed/s_6_sequence_trimmed.fastq'
output = 'test_profile'

profiler = line_profiler.LineProfiler(get_pristine_seqs.run, get_pristine_seqs.get_corrupt_seq_names, get_pristine_seqs.write_pristine)
profiler.run("get_pristine_seqs.run('%s', '%s', '%s', %s)"%(psl_file, fastq_file, output, 'False'))
profiler.print_stats()

#stats = RunningStats(in_file='stats.pkl')
#profiler = line_profiler.LineProfiler(stats.quantiles)

#profiler.run("stats.quantiles([0.05, 0.5, 0.95])")
#profiler.print_stats()