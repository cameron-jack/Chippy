#!/usr/bin/env python

from optparse import make_option
import numpy as np
from cogent.util.misc import parse_command_line_parameters

from chippy.util.inline_stats import RunningStats
from chippy.parse.fastq import FastqParser, MinimalFastqParser
import time

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"
__version__ = '0.1'

try:
    import matplotlib.pyplot as plt
except ImportError:
    # case arising on a server without graphics library
    pass

def format_long(n):
    """format a long integer"""
    n = str(int(n))
    num, rem = divmod(len(n),3)
    r = []
    if rem:
        r = [n[:rem]]

    for i in range(rem, len(n), 3):
        r += [n[i: i+3]]

    return ','.join(r)

def _num_suffix(num):
    """returns suffix for expressing proportions"""
    assert num > 1
    if 2 <= num <= 3:
        suffix = 'rd'
    else:
        suffix = 'th'
    return suffix

def make_display(stats, output_file, upper_quantile, lower_quantile, sample_freq, total):
    """produces a png file"""
    print 'Computing stats'

    lower, median, upper = stats.quantiles([0.05, 0.5, 0.95])
    # plot the results
    x = np.arange(stats.length) # length of the longest sequence
    # hard-coding figure size for the time being so they can be easily
    # inserted into a google doc
    fig = plt.figure(figsize=(6,4))
    plt.plot(x, median, x, upper, 'r--', x, lower, 'r--')
    plt.xlabel('Base Number')
    plt.ylabel('Quality Scores')

    sampled = ''
    if sample_freq > 1:
        sampled = ' (sampled every %s%s seq)' % (sample_freq,
                                            _num_suffix(sample_freq))

    plt.title('Quality Scores Across All Bases%s' % sampled)

    plt.legend(('Median', 'Quantiles [%.2f - %.2f]' % (lower_quantile,
                                                       upper_quantile)),
                                                       loc='upper right',
                                                   prop=dict(size=12))

    plt.fill_between(x, upper, lower, alpha=0.1, color='c')
    plt.setp(plt.gca(), 'ylim', [0,50])
    num = format_long(total)
    plt.text(2, 45, '%s reads' % num)
    plt.savefig(output_file)

def plot_quality(input_file, upper_quantile, lower_quantile, sample_freq,
        output_file=None, limit=None):
    if output_file is None:
        output_file = 'quality.png'

    data = None
    num = 0
    start = time.time()
    for label, seq, qual in MinimalFastqParser(input_file):
        if limit is not None and num > limit:
            break

        if data is None:
            data = dict([(chr(i+66), [0]*len(seq)) for i in range(39)])

        if num % 100000 == 0:
            print 'Working on seq %s' % format_long(num)

        num += 1
        if num % sample_freq == 0:
            for i in range(len(seq)):
                data[qual[i]][i] += 1

    end = time.time()
    print 'Took %s minutes' % ((end-start)/60)
    stats = RunningStats(data)
    make_display(stats, output_file, upper_quantile, lower_quantile,
                sample_freq, num)

# 
script_info = {}
descr = "Get Quality Information from a fastq file and plot the quality "\
        "scores across all bases."
script_info['brief_description']= descr
script_info['script_description'] = descr
script_info['version'] = '0.1.alpha'
script_info['script_usage']=[]
script_info['script_usage'].append(
    ("Example 1","""General Usage:""",
    """python plot_quality.py -i <seqs.fastq> -o <qual_plot.png>"""))

script_info['help_on_no_arguments'] = True
script_info['required_options'] = [
    make_option('-i','--input_file', help='The input fastq sequence file')]
script_info['optional_options'] = [
    make_option('-o','--output_file',
                help='The file to store the quality plot to. If not '\
                'specified, plot is stored to quality.png'),
    make_option('-u', '--upper_quantile', type='float', default=0.95,
        help='upper quantile [default: %default]'),
    make_option('-l', '--lower_quantile', type='float', default=0.05,
        help='lower quantile [default: %default]'),
    make_option('-s', '--sample_freq', type='int', default=1,
        help='sample frequency [default: %default]'),
    make_option('-L', '--limit', type='int', default=0,
        help='number of sequences to limit to [default: %default]'),
        ]


def main():
    parser, opts, args = parse_command_line_parameters(**script_info)
    plot_quality(opts.input_file,
        opts.upper_quantile,
        opts.lower_quantile,
        opts.sample_freq,
        opts.output_file,
        opts.limit or None)

if __name__ == "__main__":
    main()
