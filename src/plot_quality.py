import numpy as np
import matplotlib.pyplot as plt

from light_seq import LightSeq
from inline_stats import RunningStats
from parse_fastq import FastqParser

def plot_quality(input_file, upper_quantile, lower_quantile, output_file=None):
    
    if output_file is None:
        output_file = 'quality.png'
    
    stats = RunningStats()
    for seq in FastqParser(input_file, numeric_qual=False, make_seq=LightSeq):
        stats(seq.getNumericQuality())
    
    # obtain the mean and standard deviation accross all bases
    median = stats.quantile(0.5)
    
    upper = stats.quantile(upper_quantile)#[min(mean[i]+2*sd[i], 40) for i in range(mean.shape[0])]
    lower = stats.quantile(lower_quantile)#[max(mean[i]-2*sd[i], 2) for i in range(mean.shape[0])]
    
    # plot the results
    x = np.arange(stats.length) # length of the longest sequence
    plt.plot(x, median, x, upper, 'r--', x, lower, 'r--')
    plt.xlabel('Base Number')
    plt.ylabel('Quality Scores')
    plt.title('Quality Scores Across All Bases')
    plt.legend(('Median', 'Quantiles [%.2f - %.2f]' % (lower_quantile, upper_quantile)), loc='upper right')
    plt.fill_between(x, upper, lower, alpha=0.1, color='c')
    plt.setp(plt.gca(), 'ylim', [0,50])
    plt.savefig(output_file)

if __name__ == "__main__":
    from cogent.util.misc import parse_command_line_parameters
    from optparse import make_option
    
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
            help='lower quantile [default: %default]')]
    parser, opts, args = parse_command_line_parameters(**script_info)
    plot_quality(opts.input_file,
        opts.upper_quantile,
        opts.lower_quantile,
        opts.output_file)