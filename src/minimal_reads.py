#!/usr/bin/env python
import os
import time
import warnings
warnings.filterwarnings('ignore', 'Not using MPI as mpi4py not found')

from numpy import savez, int32, int8, uint8, array, Inf

from cogent import LoadTable
from cogent.parse.table import ConvertFields
from cogent.util.progress_display import display_wrap

from cogent.util.misc import parse_command_line_parameters
from optparse import make_option

from cogent.parse.bowtie import BowtieOutputParser
from util import create_path
from definition import NULL_STRAND, PLUS_STRAND, MINUS_STRAND

chroms = tuple(['Do All'] + map(str, range(1,20)+['X', 'Y']))

def mapped_coords(mapfile, chrom_name, limit, dry_run):
    """From a bowtie output file and a chromosome name, generate
    a list of coordinates [(start, end, strand)]."""
    
    parser = BowtieOutputParser(mapfile, row_converter=None)
    header = parser.next()
    
    all_coords = {}
    count_records = 0
    
    if dry_run:
        limit = 100
    
    for record in parser:
        chrom = record[2]
        if chrom != chrom_name:
            continue
        
        count_records += 1
        
        strand = [MINUS_STRAND, PLUS_STRAND][record[1] == '+']
        # coordinates for the match
        start = record[3]
        length = len(record[4])
        try:
            all_coords[(start, length, strand)] += 1
        except KeyError:
            all_coords[(start, length, strand)] = 1
        
        if count_records >= limit:
            break
        
    if dry_run:
        print 'Total mapped reads for %s: %d' % (chrom_name, count_records)
    all_coords = array([key + (val,) for key, val in all_coords.items()])
    all_coords = all_coords.astype(int32)
    table = LoadTable(header=['start', 'length', 'strand', 'freq'],
                        rows=all_coords)
    return table

def what_chromosomes(chrom_name, chroms=chroms):
    """returns list of chromosomes to be done"""
    chroms = list(chroms[:])
    if chrom_name == chroms[0]:
        chroms = chroms[1:]
    elif chrom_name in chrom_names:
        chroms = [strchrom_name]
    else:
        raise RuntimeError('Unknown chrom_name: %s' % str(chrom_name))
    
    return chroms


@display_wrap
def run(input_file, outdir, chroms, limit, run_record, dry_run, ui=None):
    print 'Starting'
    chroms = what_chromosomes(chroms)
    for chrom in ui.series(chroms, noun='Reading bowtie map data for chrom'):
        chrom = 'chr%s' % chrom
        table = mapped_coords(input_file, chrom, limit, dry_run)
        file_name = os.path.join(outdir, '%s.txt.gz' % chrom)
        table = table.sorted()
        end = time.time()
        if run_record:
            run_record.addMessage('minimal_reads', 'stdout',
                'Unique reads on %s' % chrom, table.Shape[0])
        
        if not dry_run:
            create_path(outdir)
            table.writeToFile(file_name, sep='\t')
        else:
            print 'will create outdir=%s' % outdir
            print 'will create outfile=%s' % file_name
        
        del(table)
    return run_record

script_info = {}
descr = "Create a coordinates tabe delimited text file for a specified"\
        " chromosome where the coordinates are all those on the chromosome."\
        "  The file contains the start, length, strand, freq of the match and is"\
        " sorted by all columns in descending order."

script_info['brief_description']= descr
script_info['script_description'] = descr
script_info['version'] = '0.1.alpha'

script_info['required_options'] = [
    make_option('-i','--input_file',
                help='The input map file (output from bowtie)'),
    make_option('-o','--outdir',
                help='Location to save the output counts.'),
    make_option('-c', '--chrom_name', 
        type='choice', default=chroms[0],
        choices= chroms,
        help="the chromosome name")
    ]

script_info['optional_options'] = [
 make_option('-d',
             '--dry_run',
             action='store_true',
             help='Whether to write data',
             default=False),
 make_option('-L', '--limit', type='int', default=Inf,
    help='number of records to read in (defaults to all)')
]


def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    
    if opts.limit == '' or opts.limit is None:
        opts.limit = Inf
    
    run(opts.input_file, opts.outdir, opts.chrom_name, opts.limit, None,
        opts.dry_run)

if __name__ == "__main__":
    main()

