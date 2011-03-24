#!/usr/bin/env python

import warnings
warnings.filterwarnings('ignore', "Not using MPI as mpi4py not found")

import numpy
from cogent import LoadTable
from cogent.util.progress_display import display_wrap

from chippy.core.read_count import CacheLaneCounts
from gene_data import GetGeneIndexes
from chippy.draw.util import smooth

try:
    from matplotlib import pyplot
    from matplotlib.ticker import MultipleLocator
except ImportError:
    # case arising on a server without graphics library
    pass

get_gene_index = GetGeneIndexes()

def get_gene_coords(infile_name, chrom_name):
    """returns TSS, strand for the nominated chromosome"""
    table = LoadTable(infile_name, sep='\t')
    table = table.filtered(lambda x: x == chrom_name, columns='CoordName')
    table = table.sorted(columns=['CoordName', 'Start'])
    rows = []
    for row in table:
        strand = row['Strand']
        if strand == -1:
            TSS = row['End']
        else:
            TSS = row['Start']

        rows += [(TSS, strand)]

    return rows

@display_wrap
def run(stable_id_file, count_dir, control_lane, treatment_lane, plot_title, ui):
    genes = map(get_gene_index,
            [l.strip() for l in open(stable_id_file).readlines()])
    num_genes = len(genes)
    treatment_cache = CacheLaneCounts(treatment_lane, count_dir)
    control_cache = CacheLaneCounts(control_lane, count_dir)
    treatments = {}
    controls = {}
    trt_binned = None
    ctl_binned = None
    for gene in genes:
        if gene.chrom in treatments:
            treatment = treatments[gene.chrom]
            control = controls[gene.chrom]
        else:
            treatment = treatment_cache[gene.chrom]
            control = control_cache[gene.chrom]
            treatments[gene.chrom] = treatment
            controls[gene.chrom] = control
            if trt_binned is None:
                trt_binned = numpy.zeros(control.shape[1], float)
                ctl_binned = numpy.zeros(control.shape[1], float)
        ctl_binned += control[gene.index]
        trt_binned += treatment[gene.index]

    trt_binned /= num_genes
    ctl_binned /= num_genes

    diff = trt_binned - ctl_binned

    sd = numpy.sqrt(trt_binned + ctl_binned)
    norm_diff = diff / sd
    norm_diff /= num_genes
    norm_diff = smooth(norm_diff, 75, 'flat')
    x = numpy.arange(-(control.shape[1]/2), (control.shape[1]/2))
    fig = pyplot.figure(figsize=(10,5))
    minor_locator = MultipleLocator(100)
    pyplot.plot(x, norm_diff)
    pyplot.xlabel('Position relative to TSS')
    pyplot.ylabel('Normalised difference (averaged by gene)')
    pyplot.title("%s: N=%s" % (plot_title, len(genes)))
    pyplot.gca().xaxis.set_minor_locator(minor_locator)
    pyplot.setp(pyplot.gca(), ylim=(-0.005,0.02))
    pyplot.grid(True)
    pyplot.savefig('%s.pdf' % stable_id_file)
    print '\n\nDone!'

if __name__ == "__main__":

    from cogent.util.misc import parse_command_line_parameters
    from optparse import make_option

    script_info = {}
    descr = "Given an R-output file with the genes that we are interested in, "\
            "a collective plot of the counts is plotted and saved. User "\
            "needs to specify the location of the counts file, the lanes for "\
            "control and treatement and location of the file which holds the "\
            "gene information."

    script_info['brief_description']= descr
    script_info['script_description'] = descr
    script_info['version'] = '0.1.alpha'
    script_info['script_usage']=[]
    script_info['script_usage'].append(
        ("Example 1","""Control lane 7; Treatment Lane 8""",
        """python segment_count.py -s stable_id_file -i count_dir -c 7 -t 8 """\
        """-p 'G2 > G1'"""))

    script_info['help_on_no_arguments'] = True
    script_info['required_options'] = [
        make_option('-s','--stable_id_file',
                    help='The R generated output file containing gene info.'),
        make_option('-i','--count_dir',
                    help='Directory where the count files are stored for '\
                    'both treatment and control'),
        make_option('-c', '--control_lane',
                    help="The lane number for the control run"),
        make_option('-t', '--treatment_lane',
                    help="The lane number for the treatment run"),
        make_option('-p', '--plot_title', help="Plot Title")
        ]

    parser, opts, args = parse_command_line_parameters(**script_info)

    run(opts.stable_id_file, opts.count_dir, opts.control_lane,
        opts.treatment_lane, opts.plot_title)
