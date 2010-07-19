import sys
try:
    import line_profiler
except ImportError:
    print 'Requires line_profiler module'
    print '$ sudo easy_install -U line_profiler'
    exit()

from parse_fastq import FastqParser

profiler = line_profiler.LineProfiler(FastqParser)
parser = FastqParser('../data/sample.txt')
profiler.run("[(n,s,q) for n, s, q in parser]")
profiler.print_stats()