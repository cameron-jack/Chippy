#!/usr/bin/env python

import os
from cogent import LoadTable

from segment_count import get_gene_coords
from make_counts import get_read_counts_bowtie, get_read_counts_sam, \
     get_file_length

def run(input_file, outdir, chrom_name, window_size, sam_output, control):

    try:
        chrom_name = int(chrom_name)
    except ValueError:
        pass

    fname_prefix = os.path.split(input_file)[1].split('.')[:-1]
    fname_prefix = '.'.join(fname_prefix)

    # static information files
    chrom_lengths_file = '../data/mouse_chrom_lengths_release_58.txt'
    gene_coords_file = '../data/mouse_gene_coords.txt'

    # read in the chromosome lengths and create a dictionary
    chrom_lengths = LoadTable(chrom_lengths_file, sep='\t')
    chrom_lengths = dict(chrom_lengths.getRawData(['chrom', 'length']))
    chrom_str = 'chr%s' % chrom_name
    mouse_gene_coords = get_gene_coords(gene_coords_file, chrom_name)
    print '\n\nGetting Counts for ' + chrom_str

    if sam_output:
        counter = get_read_counts_sam(input_file, chrom_str)
    else:
        counter = get_read_counts_bowtie(input_file, chrom_str,
                                         chrom_lengths[chrom_name])
    outfile_name = os.path.join(outdir,
            '%s-window_%s-%s.npy' % (fname_prefix, window_size, chrom_str))
    print 'Saving binary for %s to %s' % (chrom_str, outfile_name)
    counter.save(outfile_name, mouse_gene_coords,
                 window_size, control)

if __name__ == "__main__":
    from cogent.util.misc import parse_command_line_parameters
    from optparse import make_option

    script_info = {}
    descr = "Save information for each chromosome " \
            "about the number of hits per coordinate (where coordinates of " \
            "interst are the transcription start sites for that chromosome " \
            "plus/minus a set number of bases) as found by the bowtie "\
            "alignment to the reference genome, in the form of numpy arrays."

    script_info['brief_description']= descr
    script_info['script_description'] = descr
    script_info['version'] = '0.1.alpha'
    script_info['script_usage']=[]
    script_info['script_usage'].append(
        ("Example 1","""Bowtie Output with window size 2000 bps:""",
        """python region_analysis.py -i somefile.map -o outdir"""))
    script_info['script_usage'].append(
        ("Example 2","""SAM Output with window size 2000 bps:""",
        """python region_analysis.py -i somefile.map -o outdir -s"""))
    script_info['script_usage'].append(
        ("Example 3","""Bowtie Output with window size 4000 bps:""",
       """python region_analysis.py -i somefile.map -o outdir -w 4000"""))
    script_info['script_usage'].append(
        ("Example 4","""Bowtie Output for control data with window size 10000 bps:""",
       """python region_analysis.py -i somefile.map -o outdir -w 4000 -m 'control'"""))

    script_info['help_on_no_arguments'] = True
    script_info['required_options'] = [
        make_option('-i','--input_file',
                    help='The input alignment file (output from bowtie)'),
        make_option('-o','--outdir',
                    help='Location to save the output numpy arrays. The '\
                    'file basename will automatically be derived from input '\
                    'filename and appended with _chrN.npy based on the '\
                    'chromosome number in question.'),
        make_option('-c', '--chrom_name', help="the chromosome name")
        ]

    script_info['optional_options'] = [\
        make_option('-s','--sam', action='store_true',
                    dest='sam_output', default = False,
                    help='Whether the input alignment file is in SAM format '
                    +'[default: %default]'),
        make_option('-w','--window_size', type='int',
                    dest='window_size', default = 2000,
                    help='Number of bases on either site of TSS that you are '\
                    'interested in ' +'[default: %default]'),
        make_option('-m','--mode', action='store_true',
                    default=False, dest='mode',
                    help='Determines whether the data being created is for '\
                    'genomic control or otherwise')
                    ]

    parser, opts, args = parse_command_line_parameters(**script_info)

    run(opts.input_file, opts.outdir, opts.chrom_name, opts.window_size,
        opts.sam_output, opts.mode)


