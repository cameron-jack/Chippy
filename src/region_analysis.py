#!/usr/bin/env python

import os
from cogent import LoadTable

from segment_count import get_gene_coords
from make_counts import get_read_counts_bowtie, get_read_counts_sam, \
     get_file_length

def run(input_file, outfile_root, window_size, sam_output):

    # static information files
    chrom_lengths_file = '../data/mouse_chrom_lengths_release_58.txt'
    gene_coords_file = '../data/mouse_gene_coords.txt'

    # read in the chromosome lengths and create a dictionary
    chrom_lengths = LoadTable(chrom_lengths_file, sep='\t')
    chrom_lengths = dict(chrom_lengths.getRawData(['chrom', 'length']))

    for chrom in chrom_lengths:
        chrom_str = 'chr' + str(chrom)
        mouse_gene_coords = get_gene_coords(gene_coords_file, chrom)
        print '\n\nGetting Counts for ' + chrom_str

        if sam_output:
            counter = get_read_counts_sam(input_file, chrom_str)
        else:
            counter = get_read_counts_bowtie(input_file, chrom_str,
                                             chrom_lengths[chrom])
        print 'Saving binary for ' + chrom_str
        counter.save(outfile_root + '_%s'%chrom_str, mouse_gene_coords,
                     window_size)

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
        """python region_analysis.py -i somefile.map -o outfile_root"""))
    script_info['script_usage'].append(
        ("Example 2","""SAM Output with window size 2000 bps:""",
        """python region_analysis.py -i somefile.map -o outfile_root -s"""))
    script_info['script_usage'].append(
        ("Example 3","""Bowtie Output with window size 4000 bps:""",
       """python region_analysis.py -i somefile.map -o outfile_root -w 4000"""))

    script_info['help_on_no_arguments'] = True
    script_info['required_options'] = [
        make_option('-i','--input_file',
                    help='The input alignment file (output from bowtie)'),
        make_option('-o','--outfile_root',
                    help='Location and root of the output numpy arrays. The '\
                    'file basename will automatically be appended with '\
                    '_chrN.npy based on the chromosome number in question.'),
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
                    ]

    parser, opts, args = parse_command_line_parameters(**script_info)

    run(opts.input_file, opts.outfile_root, opts.window_size, opts.sam_output)


