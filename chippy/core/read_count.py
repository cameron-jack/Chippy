from __future__ import division
import sys
sys.path.extend(['..'])

import warnings
warnings.filterwarnings('ignore', 'Not using MPI as mpi4py not found')

import re
import gc
from cogent.util.progress_display import display_wrap

from chippy.util.run_record import RunRecord
from chippy.util.util import run_command
from gzip import GzipFile
from math import ceil
import numpy

__author__ = 'Cameron Jack, Gavin Huttley'
__copyright__ = 'Copyright 2011-2013, Gavin Huttley, Cameron Jack, Anuj Pahwa'
__credits__ = ['Cameron Jack', 'Gavin Huttley']
__license__ = 'GPL'
__maintainer__ = 'Cameron Jack'
__email__ = 'cameron.jack@anu.edu.au'
__status__ = 'pre-release'
__version__ = '0.2'

def get_roi_scores_from_chrom(chrom_array, chrom, all_rois, roi_ids_by_chrom):
    """ for genes in the current chrom, slice chrom_array and sum the scores """
    try:
        id_list = roi_ids_by_chrom[chrom]
    except KeyError:
        id_list = []
    for id in id_list:
        all_rois[id].score = chrom_array[all_rois[id].start:all_rois[id].end]

@display_wrap
def read_wiggle(wiggle_path, ROIs, chr_prefix='', chrom_size=300000000,
        ui=None):
    """
        Wiggles are horrible to read from as they have either fixed span
        scores per line or variable position and score per line. It's hard
        to 'bundle up' many lines before considerable CPU needs to be done
        to tie line against region.

        There is NO way to estimate the number of read tags that were used
        in the original study. We will make bases = total_score * span.
        Tags = bases/75.
    """
    rr = RunRecord('read_wiggle')
    total_score = 0 # used to calculate normalisation factors
    roi_ids_by_chrom = {} # to save time in parsing positions
    all_rois = {} # ROIs by uniqueID - converted back to list at end
    for roi in ROIs:
        if not roi.chrom in roi_ids_by_chrom.keys():
            roi_ids_by_chrom[roi.chrom] = []
        roi_ids_by_chrom[roi.chrom].append(roi.uniqueID())
        all_rois[roi.uniqueID()] = roi

    if wiggle_path.endswith('.gz'):
        wig_file = GzipFile(wiggle_path, 'rb')
    else:
        try:
            wig_file = open(wiggle_path, 'r')
        except IOError:
            rr.dieOnCritical('Could not open file', wiggle_path)

    # get total lines in wig for pacing the progress bar
    if not wiggle_path.endswith('.gz'):
        command = 'wc -l ' + wiggle_path
        returncode, stdout, stderr = run_command(command)
        if returncode:
            rr.addWarning('could not run wc to count WIG lines', 'error')
            total_lines = 1
        else:
            total_lines = int(stdout.strip().split(' ')[0])
            rr.addInfo('total lines in '+wiggle_path, total_lines)

    # Read each piece of the file into an artificial chromosome (Numpy array)
    # and slice out the gene regions that we have for each gene in that chrom

    chrom_array = numpy.zeros(chrom_size, dtype=numpy.float32)

    current_chrom = None
    for i, line in enumerate(wig_file):
        if i % 100 == 0:
            msg = 'Reading wiggle entries [' + str(i) +\
                  ' / ' + str(total_lines) + ']'
            progress = (float(i)/float(total_lines))
            ui.display(msg=msg, progress=progress)

        if line.startswith('track'):
            continue
        elif line.startswith('fixed'):
            # fixedStep chrom=chr10 start=56001 step=20 span=20
            step_type = 'fixed'
            step_parts = line.split(' ')
            step = [val.strip('step=').strip()\
                    for val in step_parts if val.startswith('step')][0]
            span = [val.strip('span=').strip()\
                    for val in step_parts if val.startswith('span')][0]
            chrom = [val.strip('chrom='+chr_prefix).strip()\
                     for val in step_parts if val.startswith('chrom')][0]

            if chrom == 'M':
                chrom = 'MT'

            if current_chrom is None:
                current_chrom = chrom
            elif current_chrom != chrom: # Empty chrom_array into genes
                get_roi_scores_from_chrom(chrom_array, chrom, all_rois,
                        roi_ids_by_chrom)
                current_chrom = chrom
                total_score += numpy.sum(chrom_array)
                chrom_array[:] = 0

            start = [val.strip('start=').strip()\
                     for val in step_parts if val.startswith('start')][0]
            pos = int(start)
            step = int(step)
            span = int(span)
        elif line.startswith('variable'):
            step_type = 'variable'
            step_parts = line.split(' ')
            chrom = [val.strip('chrom='+chr_prefix).strip()\
                     for val in step_parts if val.startswith('chrom')][0]

            if chrom == 'M':
                chrom = 'MT'

            if current_chrom is None:
                current_chrom = chrom
            elif current_chrom != chrom: # Empty chrom_array into genes
                get_roi_scores_from_chrom(chrom_array, chrom, all_rois,
                        roi_ids_by_chrom)
                current_chrom = chrom
                total_score += numpy.sum(chrom_array)
                chrom_array[:] = 0
        else:
            if step_type == 'fixed':
                chrom_array[pos:pos+span] = float(line.strip())/span
                pos += step
            else: #step_type == 'variable'
                if '\t' in line:
                    line_parts = line.split('\t')
                else:
                    line_parts = line.split(' ')
                chrom_array[int(line_parts[0])] = float(line_parts[1].strip())

    # empty chrom_array into genes_score from the final section
    get_roi_scores_from_chrom(chrom_array, chrom, all_rois,
            roi_ids_by_chrom)
    total_score += numpy.sum(chrom_array)

    ROIs = [all_rois[key] for key in all_rois.keys()]
    num_tags = total_score / 75
    num_bases = total_score
    mapped_tags = num_tags

    return ROIs, num_tags, num_bases, mapped_tags

@display_wrap
def read_BEDgraph(bedgraph_path, ROIs, chr_prefix='', ui=None):
    """
        Map BEDgraph entries to ROIs. This is essential as data uploaded
        to GEO will be BEDgraphs - we need to be able to read our own data!

        BEDgraph positions are 0-relative. We will need to add 1 to the
        last position to make it the same as Python's 0-1 system.
        Since there is no way to know the actual number of tags used in
        creating a BEDgraph we will simply assume 75-bp reads that are
        fully mapped. This will likely be conservative for most data
        currently available.
    """
    rr = RunRecord('read_BEDgraph')
    num_tags = 0; num_bases = 0; # num_tags = num_bases/75
    total_bases = 0 # for estimating total experiment tags

    try:
        if '.gz' in bedgraph_path:
            bed_graph = GzipFile(bedgraph_path, 'rb')
        else:
            bed_graph = open(bedgraph_path, 'r')
    except IOError:
        rr.dieOnCritical('No bedgraph file found', bed_graph)

    # get total lines in file for pacing the progress bar
    command = 'wc -l ' + bedgraph_path
    returncode, stdout, stderr = run_command(command)
    if returncode:
        rr.addWarning('could not run wc to count BED lines', 'error')
        total_BEDgraph_lines = 1
    else:
        total_BEDgraph_lines = int(stdout.strip().split(' ')[0])
        rr.addInfo('total lines in '+bedgraph_path, total_BEDgraph_lines)

    # separate ROIs by chrom - performance optimisation
    roi_chrom_dict = {}
    for roi in ROIs:
        if roi.chrom in roi_chrom_dict.keys():
            roi_chrom_dict[roi.chrom].append(roi)
        else:
            roi_chrom_dict[roi.chrom] = [roi]
    # sort by roi.start
    for chrom_key in roi_chrom_dict.keys():
        roi_chrom_dict[chrom_key] = sorted(roi_chrom_dict[chrom_key],
            key=lambda roi: roi.start)

    filled_ROIs = []
    for i, bed_entry in enumerate(bed_graph):
        if i % 100 == 0:
            msg = 'Reading BEDgraph entries [' + str(i) +\
                  ' / ' + str(total_BEDgraph_lines) + ']'
            progress = (float(i)/float(total_BEDgraph_lines))
            ui.display(msg=msg, progress=progress)

        if bed_entry.lower().startswith('track'):
            continue

        if bed_entry.startswith(chr_prefix):
            entry_parts = bed_entry.split('\t')
            if len(entry_parts) == 4:
                bed_chrom = str(entry_parts[0].lstrip(chr_prefix))
                try:
                    bed_start = int(entry_parts[1].strip())
                    # 'Stop' conversion to Python space, add 1.
                    bed_stop = int(entry_parts[2].strip()) + 1
                    bed_score = int(entry_parts[3].strip())
                except ValueError:
                    # Could be a browser track line
                    continue
        else:
            continue

        if not bed_chrom in roi_chrom_dict.keys():
            continue

        total_bases += bed_stop - bed_start
        for roi in roi_chrom_dict[bed_chrom]:
            if bed_stop >= roi.start: # potential for overlap
                if bed_start <= roi.end:
                    #add count to ROI
                    try:
                        roi.counts, counted_bases = roi.add_counts_to_ROI(
                                bed_start, bed_stop, bed_score)
                    except RuntimeError:
                        # input read in wrong direction or 0 sized
                        continue

                    # num_tags -> will be 'counted_bases' / 100.
                    num_bases += counted_bases
                else: # no more entries for ROI
                    filled_ROIs.append(roi)
                    # remove ROI from further consideration
                    del roi_chrom_dict[bed_chrom][0]
        # get remaining ROIs not yet in filled_ROIs
    remaining_ROIs = []
    for chrom_key in roi_chrom_dict.keys():
        for roi in roi_chrom_dict[chrom_key]:
            remaining_ROIs.append(roi)

    rr.addInfo('Filled Regions of Interest', len(filled_ROIs))
    rr.addInfo('Unfilled Regions of Interest', len(remaining_ROIs))

    # We estimate that each read is 75 after trimming.
    num_tags = int(ceil(num_bases / 75))
    # Estimated total tags is total_bases / 75
    total_tags = total_bases / 75

    return filled_ROIs + remaining_ROIs, num_tags, num_bases, \
            total_tags

@display_wrap
def read_BED(bedfile_path, ROIs, chr_prefix='', ui=None):
    """ read BED entries and add into each ROI as appropriate

    BED entries are 0-offset for start, 1-offset for end - same as Python.

    Output ROIs as well as total tags read and total bases added as
    these can be used by later stage for normalisation.
    """
    rr = RunRecord('read_BED')
    num_tags = 0; num_bases = 0

    try:
        if '.gz' in bedfile_path:
            bed_data = GzipFile(bedfile_path, 'rb')
        else:
            bed_data = open(bedfile_path, 'r')
    except IOError:
        rr.dieOnCritical('No BED file found', bedfile_path)

    # get total lines in file for pacing the progress bar, but
    # this is also the total number of mapped reads for normalisation
    command = 'wc -l ' + bedfile_path
    returncode, stdout, stderr = run_command(command)
    if returncode:
        rr.addWarning('could not run wc to count BED lines', 'error')
        total_BED_lines = 1
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
        bed_chrom = str(bed_parts[0]).lstrip(chr_prefix)
        try:
            bed_start = int(bed_parts[1].strip())
            bed_end = int(bed_parts[2].strip())
        except ValueError:
            rr.addWarning('BED file improperly formatted line', bed_entry)

        if not bed_chrom in roi_chrom_dict.keys():
            continue

        for roi in roi_chrom_dict[bed_chrom]:
            if bed_end >= roi.start: # potential for overlap
                if bed_start <= roi.end:
                    #add count to ROI
                    try:
                        roi.counts, counted_bases = roi.add_counts_to_ROI(
                                bed_start, bed_end)
                    except RuntimeError:
                        # input read in wrong direction or 0 sized
                        continue

                    num_tags += 1
                    num_bases += counted_bases
                else: # no more entries for ROI
                    filled_ROIs.append(roi)
                    # remove ROI from further consideration
                    del roi_chrom_dict[bed_chrom][0]

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
                    try:
                        roi.counts, counted_bases = roi.add_counts_to_ROI(
                                entry_start, entry_end)
                    except RuntimeError:
                        # input read in wrong direction or 0 sized
                        continue

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

def get_region_counts(BAMorBED, ROIs, chr_prefix=None, chrom_size=300000000):
    """
        Direct ROIs to BAM, BEDgraph or BED file reader.
        Also can work with Wiggle files but these are very slow.
        Return ROIs, the number of read tags, total counts and mapped tags
    """

    rr = RunRecord('get_region_counts')
    if 'bam' in BAMorBED.lower():
        filled_ROIs, num_tags, num_bases, mapped_tags =\
                read_BAM(BAMorBED, ROIs, chr_prefix)
    elif 'bedgraph' in BAMorBED.lower():
        filled_ROIs, num_tags, num_bases, mapped_tags =\
                read_BEDgraph(BAMorBED, ROIs, chr_prefix)
    elif 'bed' in BAMorBED.lower():
        filled_ROIs, num_tags, num_bases, mapped_tags =\
                read_BED(BAMorBED, ROIs, chr_prefix)
    elif 'wig' in BAMorBED.lower():
        filled_ROIs, num_tags, num_bases, mapped_tags =\
                read_wiggle(BAMorBED, ROIs, chr_prefix, chrom_size)
    else:
        rr.dieOnCritical('File not recognised as BAM, BEDgraph or BED',
                BAMorBED)

    rr.addInfo('Number of read tags counted', num_tags)
    rr.addInfo('Number of total bases counted', num_bases)
    rr.addInfo('Number of mapped tags in experiment', mapped_tags)

    return filled_ROIs, num_tags, num_bases, mapped_tags
