from __future__ import division
import sys
sys.path.extend(['..'])

import warnings
warnings.filterwarnings('ignore', 'Not using MPI as mpi4py not found')
import subprocess

import re
from cogent.util.progress_display import display_wrap

from chippy.util.definition import NULL_STRAND, PLUS_STRAND, MINUS_STRAND
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

def add_counts_to_ROI(roi, entry_start, entry_end):
    """ All coordinates are in Python 0-offset space """

    if roi.strand == PLUS_STRAND:
        offset_left = entry_start - roi.start
        offset_right = entry_end - roi.start
    else:
        offset_left = roi.end - entry_end
        offset_right = roi.end - entry_start

    assert offset_left < offset_right, 'left slicing coord must be less '+\
            'than right slicing coord'

    offset_left = 0 if offset_left < 0 else offset_left
    if offset_right >= len(roi.counts):
        offset_right = len(roi.counts)

    roi.counts[offset_left:offset_right] += 1
    return roi.counts

@display_wrap
def read_BED(bedfile_path, ROIs, chr_prefix='', rr=RunRecord(), ui=None):
    """ read BED entries and add into each ROI as appropriate

    BED entries are 0-offset for start, 1-offset for end - same as Python.
    """
    if '.gz' in bedfile_path:
        bed_data = GzipFile(bedfile_path, 'rb')
    else:
        bed_data = open(bedfile_path, 'r')

    command = 'wc -l ' + bedfile_path
    returncode, stdout, stderr = run_command(command)
    if returncode != 0:
        rr.addWarning('read_BED', 'could not run wc to count BED lines',
                'error')
    else:
        total_BED_lines = int(stdout.strip().split(' ')[0])
        rr.addInfo('read_BED', 'total lines in '+bedfile_path,
                total_BED_lines)

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
                  ' / ' + str(total_BED_lines) + ']\n'
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
                    roi.counts = add_counts_to_ROI(roi, entry_start, entry_end)
                else: # no more entries for ROI
                    filled_ROIs.append(roi)
                    # remove ROI from further consideration
                    del roi_chrom_dict[entry_chrom][0]
    # get remaining ROIs not yet in filled_ROIs
    remaining_ROIs = []
    for chrom_key in roi_chrom_dict.keys():
        for roi in roi_chrom_dict[chrom_key]:
            remaining_ROIs.append(roi)
    rr.addInfo('read_BED', 'Filled Regions of Interest', len(filled_ROIs))
    rr.addInfo('read_BED', 'Unfilled Regions of Interest', len(remaining_ROIs))

    return filled_ROIs + remaining_ROIs, rr

@display_wrap
def read_BAM(bamfile_path, ROIs, chr_prefix='', rr=RunRecord(), ui=None):
    """ get sequence reads for each ROI and add the counts in.
    Relies on Samtools view.

    Cogent entries are 0-offset with ends+1 so they slice perfectly into arrays
    SAM entries are 1-offset so substract from _start to get proper slice.
    """
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
                    rr.addWarning('read_BAM',
                            'Possibly incorrect chromosome prefix', part)
            if not reported:
                rr.addWarning('read_BAM',
                        'Possibly incorrect chromosome prefix', stderr)

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
                    entry_start = int(entry_parts[3]) - 1 # translate into 0-based space
                    length_components = map(int, re.findall(r"([\d]+)M", entry_parts[5]))
                    entry_length = sum(length_components)
                    entry_end = entry_start + entry_length
                    roi.counts = add_counts_to_ROI(roi, entry_start, entry_end)
                else:
                    invalid_flag_count += 1
            filled_ROIs.append(roi)
        else:
            rr.display()
            raise RuntimeError('samtools view failed')
    rr.addInfo('read_BAM', 'Number of BAM records evaluated', bam_lines_seen)
    rr.addInfo('read_BAM', 'Number of valid BAM records seen', valid_flag_count)
    rr.addInfo('read_BAM', 'Number of invalid BAM records seen', invalid_flag_count)
    return filled_ROIs, rr

def get_region_counts(BAMorBED, ROIs, chr_prefix=None, rr=RunRecord()):
    """ direct ROIs to BAM or BED file reader """

    if 'bam' in BAMorBED.lower():
        filled_ROIs, rr = read_BAM(BAMorBED, ROIs, chr_prefix, rr=rr)
    elif 'bed' in BAMorBED.lower():
        filled_ROIs, rr = read_BED(BAMorBED, ROIs, chr_prefix, rr=rr)
    else:
        rr.display()
        raise RuntimeError("File name given doesn't resemble BAM or BED ",
                BAMorBED)
    return filled_ROIs, rr

