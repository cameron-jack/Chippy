#!/usr/bin/env python
import os
import re
import time
import warnings
import math
warnings.filterwarnings('ignore', 'Not using MPI as mpi4py not found')

from numpy import int32, int8, uint8, array, Inf, empty

from cogent import LoadTable
from cogent.util.progress_display import display_wrap

from cogent.util.misc import parse_command_line_parameters
from optparse import make_option

from chippy.parse.sam import MinimalSamParser
from chippy.util.util import create_path
from chippy.ref.util import chroms
from chippy.util.definition import NULL_STRAND, PLUS_STRAND, MINUS_STRAND, \
    LOG_DEBUG, LOG_INFO, LOG_WARNING, LOG_ERROR, LOG_CRITICAL

__author__ = "Gavin Huttley, Cameron Jack"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley, Cameron Jack"
__credits__ = ["Gavin Huttley, Cameron Jack"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"
__version__ = '0.1'


# TODO the following is mouse specific, needs to be generalised
chroms = ('Do All',) + chroms['mouse']

def mapped_coords(samfile, min_quality, limit, dry_run):
    """From a SAM file and a chromosome name, generate
    a list of coordinates [(start, end, strand)]."""
    
    parser = MinimalSamParser(samfile)

    header = parser.next()
    all_coords = None
    aligner_chr_prefix = None
    count_records = 0
    count_skipped = 0
    
    if dry_run:
        limit = 100

    unique_str = 'XT:A:U' # 'unique' mapping identifier

    for record in parser:

        mapQ = record[4] # 5th SAM field gives mapping quality in Phred units
        optional_field = record[11] # 12th SAM field gives optional strings

        if optional_field == unique_str and min_quality <= mapQ < 255:
            count_records += 1
            chrom = record[2] # 3rd SAM field gives chromosome number
                    
            if aligner_chr_prefix is None:
                aligner_chr_prefix = re.findall(r'^[^0-9xyXY]+', chrom)[0]
                all_coords = dict(('%s%s' % (aligner_chr_prefix, chr), {})
                                for chr in chroms)

            strand = record[1]
            
            # coordinates for the match
            start = record[3] # 4th SAM field is the starting position on the chromosome
            length = len(record[9]) # 10th SAM field is the sequence
            try:
                all_coords[chrom][(start, length, strand)] += 1
            except KeyError:
                all_coords[chrom][(start, length, strand)] = 1
        else:
            count_skipped += 1

        if count_records + count_skipped >= limit:
            break
    
    # make a dict with consistent chrom keys
    all_consistent_keys = {}
    for chrom in chroms:
        aligner_chr_name = '%s%s' % (aligner_chr_prefix, chrom)
        chrom_name = 'chr%s' % chrom
        all_consistent_keys[chrom_name] = all_coords[aligner_chr_name]
    
    if dry_run:
        print 'Total mapped reads: %d' % (count_records)
        print 'Total skipped reads: %i' % (count_skipped)

    return all_consistent_keys

def make_chrom_coord_table(all_coords, chrom):
    """returns a Table of sorted coordinates for just the indicated chromosome"""
    header = ['start', 'length', 'strand', 'freq']
    coords = all_coords[chrom]
    coords = array([key + (val,) for key, val in coords.items()])
    coords = coords.astype(int32)
    if coords.shape[0] == 0:
        coords = empty((0,len(header)))
    chrom_table = LoadTable(header=header, rows=coords)
    chrom_table = chrom_table.sorted()
    return chrom_table

def what_chromosomes(chrom_name, chroms=chroms):
    """returns list of chromosomes to be done"""
    chroms = list(chroms[:])
    if chrom_name == chroms[0]:
        chroms = chroms[1:]
    elif chrom_name in chroms:
        chroms = [chrom_name]
    else:
        raise RuntimeError('Unknown chrom_name: %s' % str(chrom_name))
    
    return chroms


@display_wrap
def run(infile_name, outdir, chroms, pval_cutoff, limit, run_record, dry_run, ui=None):
    print 'Starting SAM reduction'

    if not 0.0 <= pval_cutoff <= 1.0:
        raise RuntimeError("p-value cutoff for read mapping quality"\
                           " must be between 0 and 1: %d" % pval_cutoff)

    chroms = what_chromosomes(chroms)
    min_quality = -10*math.log10(pval_cutoff) # convert from p value to Phred score
    
    start = time.time()
    all_coords = mapped_coords(infile_name, min_quality, limit, dry_run)

    for chrom in ui.series(chroms, noun='Writing SAM data for chrom'):
        chrom = 'chr%s' % chrom
        file_name = os.path.join(outdir, '%s.txt.gz' % chrom)
        chrom_table = make_chrom_coord_table(all_coords, chrom)
        if run_record:
            run_record.addMessage('reduce', LOG_INFO,
                'Unique reads on %s' % chrom, chrom_table.Shape[0])
        
        if not dry_run:
            create_path(outdir)
            chrom_table.writeToFile(file_name, sep='\t')
        else:
            print 'will create outdir=%s' % outdir
            print 'will create outfile=%s' % file_name
        
        del(chrom_table)
    
    end = time.time()
    if run_record:
        run_record.addMessage(program_name='reduce',
                error_type=LOG_INFO, message='Time taken (mins)',
                value= ((end-start)/60.))
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
                help='number of records to read in (defaults to all)'),
 make_option('-p', '--pval_cutoff',
                type='float', default=0.001,
                help='Minimum p-value for mapping quality of reads\
                [default: %default]')
]


def main():
    # opts.input_file must be of SAM format

    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    
    if opts.limit == '' or opts.limit is None:
        opts.limit = Inf
    
    run(opts.input_file, opts.outdir, opts.chrom_name, opts.pval_cutoff,
        opts.limit, None, opts.dry_run)

if __name__ == "__main__":
    main()

