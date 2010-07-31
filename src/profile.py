import sys
try:
    import line_profiler
except ImportError:
    print 'Requires line_profiler module'
    print '$ sudo easy_install -U line_profiler'
    exit()

from plot_quality import plot_quality, make_display
from inline_stats import RunningStats
from light_seq import LightSeq

stats = RunningStats(in_file='stats.pkl')
profiler = line_profiler.LineProfiler(stats.quantiles)

profiler.run("stats.quantiles([0.05, 0.5, 0.95])")
profiler.print_stats()