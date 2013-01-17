from __future__ import division
import sys
sys.path.extend(['..'])

import warnings
warnings.filterwarnings('ignore', 'Not using MPI as mpi4py not found')
import subprocess

from cogent.util.progress_display import display_wrap

from chippy.util.definition import NULL_STRAND, PLUS_STRAND, MINUS_STRAND
from chippy.util.run_record import RunRecord
from gzip import GzipFile

__author__ = "Cameron Jack"
__copyright__ = "Copyright 2011-2012, Gavin Huttley, Cameron Jack, Anuj Pahwa"
__credits__ = ['Cameron Jack']
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
    """ entries should be 1-offset """
    if entry_start == entry_end:
        return

    entry_end += 1 # adjust to Ensembl/python slice coords

    if roi.strand == PLUS_STRAND:
        offset_left = entry_start - roi.window_start
        offset_right = entry_end - roi.window_start
    else:
        offset_left = roi.window_end - entry_end
        offset_right = roi.window_end - entry_start

    #print 'left:', offset_left, 'right:', offset_right
    assert offset_left < offset_right, 'left slicing coord must be less '+\
            'than right slicing coord'

    offset_left = 0 if offset_left < 0 else offset_left
    if offset_right >= len(roi.counts):
        offset_right = len(roi.counts)
    roi.counts[offset_left:offset_right] += 1

@display_wrap
def read_BED(bedfile_path, ROIs, rr=RunRecord(), ui=None):
    """ read BED entries and add into each ROI as appropriate

    Cogent entries are 0-offset with ends+1 so they slice perfectly into arrays
    BED entries are 0-offset so add 1 to _end to slice
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

    sorted_ROIs = sorted(ROIs, key=lambda roi: roi.window_start)
    for i, bed_entry in enumerate(bed_data):
        if i % 1000 == 0:
            ui.display('Reading BED entries [' + str(i) + ', ' + \
                    str(total_BED_lines) + ' / ' + \
                    str((i/total_BED_lines)*100) + '%]')
        bed_parts = bed_entry.split('\t')
        entry_chrom = str(bed_parts[0])
        entry_start = int(bed_parts[1])+1 # 0-offset to 1-offset
        entry_end = int(bed_parts[2]) # already 1-offset in BED
        for roi in sorted_ROIs:
            if roi.chrom.lower() != entry_chrom.lower():
                continue
            if entry_end >= roi.window_start: # potential for overlap
                if entry_start > roi.window_end: # no more entries for ROI
                    sorted_ROIs = sorted_ROIs[1:] # remove ROI
                else: #add count to slice of ROI
                    add_counts_to_ROI(roi, entry_start, entry_end)
            else:
                break # bed_entry in no ROI from here

    return ROIs, rr

@display_wrap
def read_BAM(bamfile_path, ROIs, rr=RunRecord(), ui=None):
    """ get sequence reads for each ROI and add the counts in.
    Relies on Samtools view.

    Cogent entries are 0-offset with ends+1 so they slice perfectly into arrays
    SAM entries are 1-offset so substract from _start to get proper slice.
    """
    valid_flags = set([0, 16, 83, 99, 147, 163])
    for i, roi in enumerate(ROIs):
        if i % 100 == 0:
            ui.display('Reading BAM for regions of interest [' + str(i) + \
                    ', ' + str(len(ROIs)) + ' / ' +\
                    str((i/len(ROIs))*100) + '%]')
        command = 'samtools view ' + bamfile_path + ' ' + roi.chrom +\
                ':' + str(roi.window_start+1) + '-' + str(roi.window_end)
        returncode, stdout, stderr = run_command(command)
        if returncode == 0:
            bam_lines = stdout.split('\n')
            for entry in bam_lines:
                if len(entry) == 0:
                    break
                entry_parts = entry.split('\t')
                entry_flags = int(entry_parts[1])
                if entry_flags in valid_flags:
                    entry_start = int(entry_parts[3])
                    entry_length = len(entry_parts[9])
                    # end -1 because length includes the start position
                    entry_end = entry_start + entry_length -1
                    add_counts_to_ROI(roi, entry_start, entry_end)
        else:
            rr.display()
            raise RuntimeError('samtools view failed')
    return ROIs, rr

def get_region_counts(BAMorBED, ROIs, rr=RunRecord()):
    """ direct ROIs to BAM or BED file reader """

    if 'bam' in BAMorBED.lower():
        ROIs, rr = read_BAM(BAMorBED, ROIs, rr=rr)
    elif 'bed' in BAMorBED.lower():
        ROIS, rr = read_BED(BAMorBED, ROIs, rr=rr)
    else:
        rr.display()
        raise RuntimeError("File name given doesn't resemble BAM or BED ",
                BAMorBED)
    return ROIs, rr

