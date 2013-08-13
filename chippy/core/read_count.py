from __future__ import division
import sys
sys.path.extend(['..'])

import warnings
warnings.filterwarnings('ignore', 'Not using MPI as mpi4py not found')
import subprocess

import re
from cogent.util.progress_display import display_wrap

#from chippy.util.definition import NULL_STRAND, PLUS_STRAND, MINUS_STRAND
from chippy.core.region_of_interest import ROI
from chippy.util.run_record import RunRecord
from gzip import GzipFile

__author__ = "Cameron Jack, Gavin Huttley"
__copyright__ = "Copyright 2011-2012, Gavin Huttley, Cameron Jack, Anuj Pahwa"
__credits__ = ['Cameron Jack', 'Gavin Huttley']
__license__ = "GPL"
__maintainer__ = "Cameron Jack"
__email__ = "cameron.jack@anu.edu.au"
__status__ = "development"
__version__ = '0.1'

def run_command(command):
    """executes a command"""
    PIPE = subprocess.PIPE
    r = subprocess.Popen(command, shell=True, universal_newlines=True,
        stdout=PIPE, stderr=PIPE, bufsize=-1)

    stdout, stderr = r.communicate()
    returncode = r.returncode
    return returncode, stdout, stderr

@display_wrap
def read_BED(bedfile_path, ROIs, chr_prefix='', ui=None):
    """ read BED entries and add into each ROI as appropriate

    BED entries are 0-offset for start, 1-offset for end - same as Python.

    Output ROIs as well as total tags read and total bases added as
    these can be used by later stage for normalisation.
    """
    rr = RunRecord('read_BED')
    num_tags = 0; num_bases = 0

    if '.gz' in bedfile_path:
        bed_data = GzipFile(bedfile_path, 'rb')
    else:
        bed_data = open(bedfile_path, 'r')

    # get total lines in file for pacing the progress bar, but
    # this is also the total number of mapped reads for normalisation
    command = 'wc -l ' + bedfile_path
    returncode, stdout, stderr = run_command(command)
    if returncode != 0:
        rr.addWarning('could not run wc to count BED lines', 'error')
    else:
        total_BED_lines = int(stdout.strip().split(' ')[0])
        rr.addInfo('total lines in '+bedfile_path, total_BED_lines)

    # separate ROIs by chrom - performance optimisation
    roi_chrom_dict = {}
    for roi in ROIs:
        if roi.chrom in roi_chrom_dict.keys():
            roi_chrom_dict[roi.chrom].append(roi)
        else:
            roi_chrom_dict[roi.chrom] = [roi]

    for chrom_key in roi_chrom_dict.keys():
        roi_chrom_dict[chrom_key] = sorted(roi_chrom_dict[chrom_key],
                key=lambda roi: roi.start)
    filled_ROIs = []
    for i, bed_entry in enumerate(bed_data):
        if i % 100 == 0:
            msg = 'Reading BED entries [' + str(i) +\
                  ' / ' + str(total_BED_lines) + ']'
            progress = (float(i)/float(total_BED_lines))
            ui.display(msg=msg, progress=progress)

        bed_parts = bed_entry.split('\t')
        entry_chrom = str(bed_parts[0]).lstrip(chr_prefix)
        entry_start = int(bed_parts[1])
        entry_end = int(bed_parts[2])

        if not entry_chrom in roi_chrom_dict.keys():
            continue

        for roi in roi_chrom_dict[entry_chrom]:
            if entry_end >= roi.start: # potential for overlap
                if entry_start <= roi.end:
                    #add count to ROI
                    roi.counts, counted_bases = roi.add_counts_to_ROI(
                            entry_start, entry_end)
                    num_tags += 1
                    num_bases += counted_bases
                else: # no more entries for ROI
                    filled_ROIs.append(roi)
                    # remove ROI from further consideration
                    del roi_chrom_dict[entry_chrom][0]
    # get remaining ROIs not yet in filled_ROIs
    remaining_ROIs = []
    for chrom_key in roi_chrom_dict.keys():
        for roi in roi_chrom_dict[chrom_key]:
            remaining_ROIs.append(roi)
    rr.addInfo('Filled Regions of Interest', len(filled_ROIs))
    rr.addInfo('Unfilled Regions of Interest', len(remaining_ROIs))

    return filled_ROIs + remaining_ROIs, num_tags, num_bases, total_BED_lines

@display_wrap
def read_BAM(bamfile_path, ROIs, chr_prefix='', ui=None):
    """ get sequence reads for each ROI and add the counts in.
    Relies on Samtools view.

    Cogent entries are 0-offset with ends+1 so they slice perfectly into arrays
    SAM entries are 1-offset so substract from _start to get proper slice.

    Output ROIs as well as total tags read and total bases added as
    these can be used by later stage for normalisation.
    """
    rr = RunRecord('read_BAM')
    # Collect stats for counts normalisation
    num_tags = 0
    num_bases = 0
    mapped_tags = 0

    # get total number of mapped tags first
    command = 'samtools idxstats ' + bamfile_path
    returncode, stdout, stderr = run_command(command)
    if returncode != 0:
        rr.dieOnCritical('Samtools idxstats died', 'Indexed correctly?')
    else:
        lines = stdout.split('\n')
        for line in lines:
            line_parts = line.split('\t')
            if len(line_parts) == 4:
                mapped_reads = line_parts[2]
                mapped_tags += int(mapped_reads)

    valid_flags = set([0, 16, 83, 99, 147, 163])
    filled_ROIs = []
    bam_lines_seen = 0
    valid_flag_count = 0
    invalid_flag_count = 0
    for i, roi in enumerate(ui.series(ROIs, 'Loading Regions of Interest')):
        chrom = chr_prefix + roi.chrom # default is ''
        command = 'samtools view ' + bamfile_path + ' ' + chrom +\
                ':' + str(roi.start+1) + '-' + str(roi.end)
        returncode, stdout, stderr = run_command(command)
        if 'fail to determine the sequence name' in stderr:
            reported = False
            for part in stderr.split(' '):
                if chr_prefix in part:
                    reported = True
                    rr.addWarning('Possibly incorrect chromosome prefix',
                            part)
            if not reported:
                rr.addWarning('Possibly incorrect chromosome prefix', stderr)

        if returncode == 0:
            bam_lines = stdout.split('\n')
            for entry in bam_lines:
                if not len(entry):
                    continue
                bam_lines_seen += 1
                entry_parts = entry.split('\t')
                entry_flag = int(entry_parts[1])
                if entry_flag in valid_flags:
                    valid_flag_count += 1
                    # translate into 0-based space
                    entry_start = int(entry_parts[3]) - 1
                    length_components = map(int, re.findall(r"([\d]+)M",
                            entry_parts[5]))
                    entry_length = sum(length_components)
                    entry_end = entry_start + entry_length
                    roi.counts, counted_bases = roi.add_counts_to_ROI(
                            entry_start, entry_end)
                    num_tags += 1
                    num_bases += counted_bases
                else:
                    invalid_flag_count += 1
            filled_ROIs.append(roi)
        else:
            rr.dieOnCritical('Samtools view failed', 'Sorted? Indexed?')

    rr.addInfo('Number of BAM records evaluated', bam_lines_seen)
    rr.addInfo('Number of valid BAM records seen', valid_flag_count)
    rr.addInfo('Number of invalid BAM records seen', invalid_flag_count)

    return filled_ROIs, num_tags, num_bases, mapped_tags

def get_region_counts(BAMorBED, ROIs, chr_prefix=None):
    """ direct ROIs to BAM or BED file reader.
    Return ROIs, the number of read tags, the total of all counts """
    rr = RunRecord('get_region_counts')
    if 'bam' in BAMorBED.lower():
        filled_ROIs, num_tags, num_bases, mapped_tags =\
                read_BAM(BAMorBED, ROIs, chr_prefix)
    elif 'bed' in BAMorBED.lower():
        filled_ROIs, num_tags, num_bases, mapped_tags =\
                read_BED(BAMorBED, ROIs, chr_prefix)
    else:
        rr.dieOnCritical('File not recognised as BAM or BED', BAMorBED)

    rr.addInfo('Number of read tags counted', num_tags)
    rr.addInfo('Number of total bases counted', num_bases)
    rr.addInfo('Number of mapped tags in experiment', mapped_tags)

    return filled_ROIs, num_tags, num_bases, mapped_tags
