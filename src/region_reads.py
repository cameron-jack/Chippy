#!/usr/bin/env python

import os
from numpy import save
from cogent import LoadTable
from cogent.parse.bowtie import BowtieOutputParser

from segment_count import get_gene_coords

from cogent.util.misc import parse_command_line_parameters
from optparse import make_option


def generate_coord_range(coords_file, window_size, chrom, chrom_length):
    """Create a list of coordinates around the TSS of the form
    [(x1, y1), (x2, y2),...]
    where yi - xi = 2*window_size"""

    tss_list = [tss for (tss,strand) in
                get_gene_coords(coords_file, chrom)]

    coords = []
    for tss in tss_list:
        start = max([tss - window_size, 0])
        end = min([tss + window_size, chrom_length])
        coords.append((start, end))
    return coords

def get_matching_coords(mapfile, chrom_str, ref_coords):
    """From a bowtie output file and a list of reference coordinates, generate
    a list of coordinates where bowtie was successfully able to map to the
    genome, as long as those coordinates are within the reference range."""

    parser = BowtieOutputParser(mapfile)
    header = parser.next()

    store_coords = []
    count_records = 0
    count_matches = 0
    for record in parser:

        chrom_name = record[2]
        if chrom_str != chrom_name:
            continue

        count_records += 1
        # save strand information as 0,1
        strand = record[1]
        if strand is '+':
            strand = 0
        else:
            strand = 1

        # coordinates for the match
        start = record[3]
        end = start+len(record[4])

        for (ref_start, ref_end) in ref_coords:
            if ref_start <= start <= ref_end and ref_start <= end <= ref_end:
                count_matches += 1
                store_coords.append([start, end, strand])
                break

    print 'Total records for %s: %d' % (chrom_str, count_records)
    print 'Total matches for %s: %d' % (chrom_str, count_matches)
    return store_coords

def run(input_file, outdir, chrom_name, window_size):
    """Given a bowtie output file in the default bowtie format obtain the
    actual reads that bowtie was able to successfully map within the coordinates
    TSS-window_size, TSS+window_size"""

    try:
        chrom_name = int(chrom_name)
    except ValueError:
        pass

    # Create a dictionary of Chromosome Lengths
    chrom_lengths_file = '../data/mouse_chrom_lengths_release_58.txt'
    chrom_lengths = LoadTable(chrom_lengths_file, sep='\t')
    chrom_lengths = dict(chrom_lengths.getRawData(['chrom', 'length']))

    # Generate the coordinate ranges
    gene_coords_file = '../data/mouse_gene_coords.txt'
    coord_range = generate_coord_range(gene_coords_file, window_size,
                                       chrom_name,chrom_lengths[chrom_name])

    chrom_str = 'chr%s' % chrom_name
    match_coords = get_matching_coords(input_file, chrom_str, coord_range)

    fn_prefix = os.path.split(input_file)[1].split('.')[0]
    fn_prefix = os.path.join(outdir, '%s_%s'%(fn_prefix, chrom_str))
    save(fn_prefix, match_coords)

# 
script_info = {}
descr = "Create a coordinates binary numpy file for a specified chromosome"\
        " where the coordinates are found in a specified window around the"\
        " TSS for genes found on that chromosome. The binary file contains"\
        " the start and stop coordinates of the match, and the strand "\
        "associated with the match"

script_info['brief_description']= descr
script_info['script_description'] = descr
script_info['version'] = '0.1.alpha'
script_info['script_usage']=[]
script_info['script_usage'].append(
    ("Example 1","""Chromosome 1 with window size 2000 bps:""",
    """python region_reads.py -i somefile.map -o outdir -c 1"""))
script_info['script_usage'].append(
    ("Example 2","""Chromosome X with window size 5000 bps:""",
    """python region_reads.py -i somefile.map -o outdir -c X -w 5000"""))

script_info['help_on_no_arguments'] = True
script_info['required_options'] = [
    make_option('-i','--input_file',
                help='The input map file (output from bowtie)'),
    make_option('-o','--outdir',
                help='Location to save the output numpy arrays. The '\
                'file basename will automatically be derived from input '\
                'filename and appended with _chrN.npy based on the '\
                'chromosome number specified.'),
    make_option('-c', '--chrom_name', help="the chromosome name")
    ]

script_info['optional_options'] = [\
    make_option('-w','--window_size', type='int',
                dest='window_size', default = 2000,
                help='Number of bases on either site of TSS that you are '\
                'interested in ' +'[default: %default]')
                ]

def main():
    """executes the script"""
    parser, opts, args = parse_command_line_parameters(**script_info)
    run(opts.input_file, opts.outdir, opts.chrom_name, opts.window_size)

if __name__ == "__main__":
    main()





