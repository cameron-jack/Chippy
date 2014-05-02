from __future__ import division
import sys
sys.path.extend(['..'])

import warnings
warnings.filterwarnings('ignore', 'Not using MPI as mpi4py not found')

import re
from cogent.util.progress_display import display_wrap

from chippy.util.run_record import RunRecord
from chippy.util.util import run_command
from gzip import GzipFile
from math import ceil
import numpy
import uuid # for generating random file names for temporary SAM entries
import os

NO_CYVCF = False
try:
    import cyvcf
except ImportError:
    cyvcf = None
    NO_CYVCF = True

__author__ = 'Cameron Jack, Gavin Huttley'
__copyright__ = 'Copyright 2011-2013, Gavin Huttley, Cameron Jack, Anuj Pahwa'
__credits__ = ['Cameron Jack', 'Gavin Huttley']
__license__ = 'GPL'
__maintainer__ = 'Cameron Jack'
__email__ = 'cameron.jack@anu.edu.au'
__status__ = 'pre-release'
__version__ = '0.2'

def get_roi_scores_from_chrom(chrom_array, chrom, all_rois, roi_ids_by_chrom):
    """ for genes in the current chrom, slice chrom_array and add to counts """
    try:
        id_list = roi_ids_by_chrom[chrom]
    except KeyError:
        id_list = []
    for id in id_list:
        if all_rois[id].strand == 1:
            all_rois[id].counts +=\
                    chrom_array[all_rois[id].start:all_rois[id].end]
        else: # negative strand genes start at the far end and come back
            start = all_rois[id].start-1
            if start < 0: # this is required to include the 0-indexed position
                all_rois[id].counts +=\
                        chrom_array[all_rois[id].end-1::-1]
            else:
                all_rois[id].counts +=\
                        chrom_array[all_rois[id].end-1:start:-1]

@display_wrap
def read_vcf(vcf_path, ROIs, chr_prefix='', chrom_size=300000000,
             ui=None):
    """
        Either use cyvcf for ultra-fast VCF reading or do it ourselves
    """
    rr = RunRecord('read_wiggle')
    total_score = 0 # used to calculate normalisation factors
    roi_ids_by_chrom = {} # to save time in parsing positions
    all_rois = {} # ROIs by uniqueID - converted back to list at end
    for roi in ROIs:
        if not roi.chrom in roi_ids_by_chrom.keys():
            roi_ids_by_chrom[roi.chrom] = []
        roi_ids_by_chrom[roi.chrom].append(roi.unique_id)
        all_rois[roi.unique_id] = roi

    if NO_CYVCF:
        # Parse the VCF ourselves
        if vcf_path.endswith('.gz'):
            vcf_file = GzipFile(vcf_path, 'rb')
        else:
            try:
                vcf_file = open(vcf_path, 'r')
            except IOError:
                rr.dieOnCritical('Could not open file', vcf_path)

        # get total lines in wig for pacing the progress bar
        if not vcf_path.endswith('.gz'):
            command = 'wc -l ' + vcf_path
        else:
            command = 'gunzip -c ' + vcf_path + ' | wc -l'

        returncode, stdout, stderr = run_command(command)
        if returncode:
            rr.addWarning('could not run wc to count WIG lines', 'error')
            total_lines = 1
        else:
            total_lines = int(stdout.strip().split(' ')[0])
            rr.addInfo('total lines in ' + vcf_path, total_lines)

    else:
        vcf_file = cyvcf.Reader(open(vcf_path, 'rb'))
        total_lines = len(vcf_file)

    # Read each piece of the file into an artificial chromosome (Numpy array)
    # and slice out the gene regions that we have for each gene in that chrom
    chrom_array = numpy.zeros(chrom_size, dtype=numpy.float32)

    current_chrom = None
    for i, line in enumerate(vcf_file):
        if i % 100 == 0:
            msg = 'Reading VCF entries [' + str(i) +\
                   ' / ' + str(total_lines) + ']'
            progress = (float(i)/float(total_lines))
            ui.display(msg=msg, progress=progress)

        if line.startswith('#'):
            continue # header line

        line_parts = line.split('\t')
        chrom = line_parts[0].strip().lstrip(chr_prefix)
        if current_chrom is None:
            current_chrom = chrom
        elif current_chrom != chrom:
            get_roi_scores_from_chrom(chrom_array, current_chrom, all_rois,
                    roi_ids_by_chrom)
            current_chrom = chrom
            chrom_array[:] = 0

        position = int(line_parts[1].strip())
        if position == 0:
            continue # telomere

        try:
            chrom_array[position:position+1] += 1.0
        except IndexError:
            rr.addWarning('position larger than expect chrom size', position)
        total_score += 1

    # Get last chromosome entries
    get_roi_scores_from_chrom(chrom_array, current_chrom, all_rois,
            roi_ids_by_chrom)

    ROIs = [all_rois[key] for key in all_rois.keys()]
    num_tags = total_score
    num_bases = total_score * 75
    mapped_tags = total_score

    return ROIs, num_tags, num_bases, mapped_tags

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
        roi_ids_by_chrom[roi.chrom].append(roi.unique_id)
        all_rois[roi.unique_id] = roi

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
    else:
        command = 'gunzip -c ' + wiggle_path + ' | wc -l'

    returncode, stdout, stderr = run_command(command)
    if returncode:
        rr.addWarning('could not run wc to count WIG lines', 'error')
        total_lines = 1
    else:
        total_lines = int(stdout.strip().split(' ')[0])
        rr.addInfo('total lines in ' + wiggle_path, total_lines)

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
                get_roi_scores_from_chrom(chrom_array, current_chrom, all_rois,
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
                get_roi_scores_from_chrom(chrom_array, current_chrom, all_rois,
                        roi_ids_by_chrom)
                current_chrom = chrom
                total_score += numpy.sum(chrom_array)
                chrom_array[:] = 0
        else:
            if step_type == 'fixed':
                try:
                    chrom_array[pos:pos+span] = float(line.strip())/span
                except IndexError:
                    rr.addWarning('position larger than expect chrom size', pos)
                pos += step
            else: #step_type == 'variable'
                if '\t' in line:
                    line_parts = line.split('\t')
                else:
                    line_parts = line.split(' ')
                pos = int(line_parts[0])
                try:
                    chrom_array[pos] = float(line_parts[1].strip())
                except IndexError:
                    rr.addWarning('position larger than expect chrom size', pos)

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
def read_BEDgraph(bedgraph_path, ROIs, chr_prefix='',
        chrom_size=300000000, ui=None):
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

    try:
        if '.gz' in bedgraph_path:
            bed_graph = GzipFile(bedgraph_path, 'rb')
        else:
            bed_graph = open(bedgraph_path, 'r')
    except IOError:
        rr.dieOnCritical('No bedgraph file found', bed_graph)

    # get total lines in file for pacing the progress bar
    if not bedgraph_path.endswith('.gz'):
        command = 'wc -l ' + bedgraph_path
    else:
        command = 'gunzip -c ' + bedgraph_path + ' | wc -l'

    returncode, stdout, stderr = run_command(command)
    if returncode:
        rr.addWarning('could not run wc to count BED lines', 'error')
        total_BEDgraph_lines = 1
    else:
        total_BEDgraph_lines = int(stdout.strip().split(' ')[0])
        rr.addInfo('total lines in '+bedgraph_path, total_BEDgraph_lines)

    roi_ids_by_chrom = {} # to save time in parsing positions
    all_rois = {} # ROIs by uniqueID - converted back to list at end
    for roi in ROIs:
        if not roi.chrom in roi_ids_by_chrom.keys():
            roi_ids_by_chrom[roi.chrom] = []
        roi_ids_by_chrom[roi.chrom].append(roi.unique_id)
        all_rois[roi.unique_id] = roi

    # Read each piece of the file into an artificial chromosome (Numpy array)
    # and slice out the gene regions that we have for each gene in that chrom
    chrom_array = numpy.zeros(chrom_size, dtype=numpy.float32)

    current_chrom = None
    for i, bed_entry in enumerate(bed_graph):
        if i % 100 == 0:
            msg = 'Reading BEDgraph entries [' + str(i) +\
                  ' / ' + str(total_BEDgraph_lines) + ']'
            progress = (float(i)/float(total_BEDgraph_lines))
            ui.display(msg=msg, progress=progress)

        if bed_entry.lower().startswith('track'):
            continue

        # check formatting with tabs first, then space
        bed_parts = bed_entry.split('\t')
        if len(bed_parts) != 4:
            bed_parts = bed_entry.split(' ')
            if len(bed_parts) != 4:
                continue # some sort of track or browser line

        bed_chrom = str(bed_parts[0]).lstrip(chr_prefix)

        # check if it's time to flush the chrom counts array
        if bed_chrom != current_chrom and current_chrom is not None:
            get_roi_scores_from_chrom(chrom_array, current_chrom,
                    all_rois, roi_ids_by_chrom)
            chrom_array[:] = 0

        current_chrom = bed_chrom
        try:
            bed_start = int(bed_parts[1].strip())
            bed_end = int(bed_parts[2].strip()) + 1 # slice offset
            # score = number of reads at this location
            bed_score = float(bed_parts[3].strip())

            try:
                chrom_array[bed_start:bed_end] += bed_score
            except IndexError:
                rr.addWarning('position larger than expect chrom size',
                        bed_start)
            num_bases += (bed_end - bed_start) * bed_score
            num_tags += ((bed_end - bed_start) * bed_score)/75
        except ValueError:
            rr.addWarning('BED file improperly formatted line', bed_entry)

    # Clear last entries
    if current_chrom is not None:
        get_roi_scores_from_chrom(chrom_array, current_chrom,
                all_rois, roi_ids_by_chrom)

    # We estimate that each read is 75 after trimming.
    num_tags = int(ceil(num_bases / 75))
    # Estimated total tags is total_bases / 75
    total_tags = num_bases / 75

    ROIs = [all_rois[key] for key in all_rois.keys()]
    return ROIs, num_tags, num_bases, total_tags

@display_wrap
def read_BED(bedfile_path, ROIs, chr_prefix='',
        chrom_size=300000000, ui=None):
    """
        Read BED entries into chromosome array. Slice into each ROI.

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
    if not bedfile_path.endswith('.gz'):
        command = 'wc -l ' + bedfile_path
    else:
        command = 'gunzip -c ' + bedfile_path + ' | wc -l'

    returncode, stdout, stderr = run_command(command)
    if returncode:
        rr.addWarning('could not run wc to count BED lines', 'error')
        total_BED_lines = 1
    else:
        total_BED_lines = int(stdout.strip().split(' ')[0])
        rr.addInfo('total lines in '+bedfile_path, total_BED_lines)

    roi_ids_by_chrom = {} # to save time in parsing positions
    all_rois = {} # ROIs by uniqueID - converted back to list at end
    for roi in ROIs:
        if not roi.chrom in roi_ids_by_chrom.keys():
            roi_ids_by_chrom[roi.chrom] = []
        roi_ids_by_chrom[roi.chrom].append(roi.unique_id)
        all_rois[roi.unique_id] = roi

    # Read each piece of the file into an artificial chromosome (Numpy array)
    # and slice out the gene regions that we have for each gene in that chrom
    chrom_array = numpy.zeros(chrom_size, dtype=numpy.float32)

    current_chrom = None
    for i, bed_entry in enumerate(bed_data):
        if i % 100 == 0:
            msg = 'Reading BED entries [' + str(i) +\
                  ' / ' + str(total_BED_lines) + ']'
            progress = (float(i)/float(total_BED_lines))
            ui.display(msg=msg, progress=progress)

        bed_parts = bed_entry.split('\t')
        if len(bed_parts) <= 3:
            continue # some sort of track or browser line

        bed_chrom = str(bed_parts[0]).lstrip(chr_prefix)

        # check if it's time to flush the chrom counts array
        if bed_chrom != current_chrom and current_chrom is not None:
            get_roi_scores_from_chrom(chrom_array, current_chrom,
                    all_rois, roi_ids_by_chrom)
            chrom_array[:] = 0

        current_chrom = bed_chrom
        try:
            bed_start = int(bed_parts[1].strip())
            bed_end = int(bed_parts[2].strip())

            # don't use the quality score
            try:
                bed_score = float(bed_parts[4].strip())
            except ValueError, TypeError:
                bed_score = 1

            try:
                chrom_array[bed_start:bed_end] += 1.0
            except IndexError:
                rr.addWarning('position larger than expect chrom size',
                        bed_start)
            num_bases += bed_end - bed_start
            num_tags += (bed_end - bed_start)/75
        except ValueError:
            rr.addWarning('BED file improperly formatted line', bed_entry)

    # Clear last entries
    if current_chrom is not None:
        get_roi_scores_from_chrom(chrom_array, current_chrom,
                all_rois, roi_ids_by_chrom)

    ROIs = [all_rois[key] for key in all_rois.keys()]
    return ROIs, num_tags, num_bases, total_BED_lines

@display_wrap
def read_BAM(bam_path, ROIs, chr_prefix='', chrom_size=300000000, ui=None):
    """
        Get sequence reads for each chrom.
        Relies on Samtools view.

        Cogent entries are 0-offset with ends+1 so they slice perfectly into arrays
        SAM entries are 1-offset so subtract from _start to get proper slice.

        Output ROIs as well as total tags read and total bases added as
        these can be used by later stage for normalisation.
    """
    rr = RunRecord('read_BAM')

    roi_ids_by_chrom = {} # to save time in parsing positions
    all_rois = {} # ROIs by uniqueID - converted back to list at end
    for roi in ROIs:
        if not roi.chrom in roi_ids_by_chrom.keys():
            roi_ids_by_chrom[roi.chrom] = []
        roi_ids_by_chrom[roi.chrom].append(roi.unique_id)
        all_rois[roi.unique_id] = roi

    # Collect stats for counts normalisation
    num_tags = 0
    num_bases = 0
    mapped_tags = 0

    # get total number of mapped tags first
    command = 'samtools idxstats ' + bam_path
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
    bam_lines_seen = 0
    valid_flag_count = 0
    invalid_flag_count = 0

    # Read each piece of the file into an artificial chromosome (Numpy array)
    # and slice out the gene regions that we have for each gene in that chrom
    chrom_array = numpy.zeros(chrom_size, dtype=numpy.float32)

    for chrom_key in roi_ids_by_chrom.keys():
        if chr_prefix not in chrom_key:
            chrom = chr_prefix + chrom_key # default is ''
        else:
            chrom = chrom_key
        chrom_fn = str(uuid.uuid4())
        # get all reads per chrom
        command = 'samtools view ' + bam_path + ' ' + chrom + ' > ' + chrom_fn
        returncode, stdout, stderr = run_command(command)

        if 'fail to determine the sequence name' in stderr:
            reported = False
            for part in stderr.split(' '):
                if chr_prefix in part:
                    reported = True
                    rr.addWarning('Possibly incorrect chromosome prefix', part)
            if not reported:
                rr.addWarning('Possibly incorrect chromosome prefix', stderr)

        if returncode == 0:
            with open(chrom_fn, 'r') as bam_lines:
                for entry in bam_lines:
                    if not len(entry):
                        continue
                    bam_lines_seen += 1

                    if bam_lines_seen % 100 == 0:
                        msg = 'Reading BAM entries [' + str(bam_lines_seen) +\
                                ' / ' + str(mapped_tags) + ']'
                        progress = (float(bam_lines_seen)/float(mapped_tags))
                        ui.display(msg=msg, progress=progress)

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

                        chrom_array[entry_start:entry_end] += 1.0

                        num_tags += 1
                        num_bases += entry_length
                    else:
                        invalid_flag_count += 1
        else:
            os.remove(chrom_fn)
            rr.dieOnCritical('Samtools view failed; Sorted? Indexed?',
                    returncode)

        os.remove(chrom_fn) # clean up temp chrom file
        get_roi_scores_from_chrom(chrom_array, chrom_key, all_rois, roi_ids_by_chrom)
        chrom_array[:] = 0

    rr.addInfo('Number of BAM records evaluated', bam_lines_seen)
    rr.addInfo('Number of valid BAM records seen', valid_flag_count)
    rr.addInfo('Number of invalid BAM records seen', invalid_flag_count)

    ROIs = [all_rois[key] for key in all_rois.keys()]
    return ROIs, num_tags, num_bases, mapped_tags

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
    elif 'vcf' in BAMorBED.lower():
        filled_ROIs, num_tags, num_bases, mapped_tags =\
                read_vcf(BAMorBED, ROIs, chr_prefix, chrom_size)
    else:
        rr.dieOnCritical('File not recognised as BAM, BEDgraph,'+\
                'BED, WIG or VCF', BAMorBED)

    rr.addInfo('Number of read tags counted', num_tags)
    rr.addInfo('Number of total bases counted', num_bases)
    rr.addInfo('Number of mapped tags in experiment', mapped_tags)

    return filled_ROIs, num_tags, num_bases, mapped_tags
