#!/usr/bin/env python
import os, sys, re, time
import subprocess

import pristine_paired

from chippy.util.definition import LOG_DEBUG, LOG_INFO, LOG_WARNING, \
    LOG_ERROR, LOG_CRITICAL

__author__ = "Gavin Huttley, Cameron Jack, Aaron Chuah"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley, Cameron Jack, Aaron Chuah"
__credits__ = ["Gavin Huttley", "Cameron Jack", "Aaron Chuah"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"
__version__ = '0.1'


def run_command(command, test=None):
    """executes a command"""
    PIPE = subprocess.PIPE
    r = subprocess.Popen(command, shell=True, universal_newlines=True,
            stdout=PIPE, stderr=PIPE, bufsize=-1)
    
    stdout, stderr = r.communicate()
    returncode = r.returncode
    return returncode, stdout, stderr

def run_blat(blat_adapters, query_file, psl_out, run_record, test):
    """run blat to remove adapter sequences"""
    command = "blat %s %s %s" % (blat_adapters, query_file, psl_out)
    if test:
        print "=== The command ==="
        print command
        return
    start = time.time()
    returncode, stdout, stderr = run_command(command, test)
    end = time.time()
    run_record.addMessage(program_name='blat',
            error_type=LOG_INFO, message='Time taken (mins)',
            value=((end-start)/60.))
    if stdout:
        print
        print ''.join(stdout)
    
    if stderr:
        print
        print ''.join(stderr)
    return run_record

def run_fastx_clipper(blat_adapters, fastq_in_fn, fastq_out_fn, run_record, num_threads, test):
    """run fastx_clipper to remove adapter sequences"""

    file_in = open (blat_adapters)
    # ignore header
    line = file_in.readline()
    # get first sequence
    adapter = file_in.readline().strip()
    file_in.close()
    decompress_threads=2
    if num_threads<8: decompress_threads=1
    compress_threads=num_threads-decompress_threads-1
    if compress_threads<1: compress_threads=1
    command = 'pigz -d -c -p %(dthreads)s %(infile)s | fastx_clipper -a %(adapter)s | pigz -p %(cthreads)s > %(outfile)s' % \
              dict(dthreads=decompress_threads, cthreads=compress_threads, infile=fastq_in_fn, adapter=adapter, 
                    outfile=fastq_out_fn)
    
    if test:
        print "=== The command ==="
        print command
        return run_record
    
    start = time.time()
    
    returncode, stdout, stderr = run_command(command, test)
    end = time.time()
    run_record.addMessage(program_name=command,
            error_type=LOG_INFO, message='Time taken (mins)',
            value=((end-start)/60.))
    if stdout:
        print
        print ''.join(stdout)

    if stderr:
        print
        print ''.join(stderr)
    return run_record

def run_fastq_qual_trim(fastq_in_fn, fastq_out_fn, run_record, num_threads, test):
    """run fastq_quality_trimmer to remove adapter sequences. Note that
       we are being reasonably conservative with the minimum length we
       pass on for alignment"""
    decompress_threads=2
    if num_threads<8: decompress_threads=1
    compress_threads=num_threads-decompress_threads-1
    if compress_threads<1: compress_threads=1
    command = 'pigz -d -c -p %(dthreads)s %(infile)s | fastq_quality_trimmer -t 3 -l 19 | pigz -p %(cthreads)s > %(outfile)s' \
              % dict(dthreads=decompress_threads, cthreads=compress_threads, infile=fastq_in_fn, outfile=fastq_out_fn)
    if test:
        print "=== The command ==="
        print command
        return run_record
    
    start = time.time()
    
    returncode, stdout, stderr = run_command(command, test)
    end = time.time()
    run_record.addMessage(program_name=command,
            error_type=LOG_INFO, message='Time taken (mins)',
            value=((end-start)/60.))
    if stdout:
        print
        print ''.join(stdout)

    if stderr:
        print
        print ''.join(stderr)
    return run_record


def run_fastx_clip_and_trim(blat_adapters, fastq_in_fn, fastq_out_fn, run_record, num_threads, test):
    """run fastx_clipper to remove adapter sequences"""

    file_in = open (blat_adapters)
    # ignore header
    line = file_in.readline()
    # get first sequence
    adapter = file_in.readline().strip()
    file_in.close()
    decompress_threads=num_threads/4
    if decompress_threads<1: decompress_threads=1
    compress_threads=num_threads-decompress_threads
    if compress_threads<1: compress_threads=1
    command = 'pigz -d -c -p %(dthreads)s %(infile)s | parallel -k -d --pipe --spreadstdin --recstart "@" -j %(threads)s "fastx_clipper -a %(adapter)s | fastq_quality_trimmer -t 3 -l 19" | pigz -p %(cthreads)s > %(outfile)s' % \
              dict(dthreads=decompress_threads, cthreads=compress_threads, threads=num_threads, infile=fastq_in_fn, adapter=adapter, outfile=fastq_out_fn)
    
    if test:
        print "=== The command ==="
        print command
        return run_record
    
    start = time.time()
    returncode, stdout, stderr = run_command(command, test)
    end = time.time()
    run_record.addMessage(program_name=command,
            error_type=LOG_INFO, message='Time taken (mins)',
            value=((end-start)/60.))
    if stdout:
        print
        print ''.join(stdout)

    if stderr:
        print
        print ''.join(stderr)
    return run_record


def run_pristine_paired(fastq_in1_fn, fastq_in2_fn, fastq_out1_fn, fastq_out2_fn, 
run_record, num_threads, test):
    """run pristine_paired.py to ensure paired filtered reads"""

    command = 'pristine_paired.py %s %s %s %s %d' % (fastq_in1_fn, fastq_in2_fn, 
fastq_out1_fn, fastq_out2_fn, num_threads)

    if test:
        print "=== The command ==="
        print command
        return run_record
    
    start = time.time()
    pristine_paired.make_pristine_pair(fastq_in1_fn, fastq_in2_fn, fastq_out1_fn, fastq_out2_fn, num_threads)
    end = time.time()
    run_record.addMessage(program_name=command,
            error_type=LOG_INFO, message='Time taken (mins)',
            value=((end-start)/60.))
    return run_record


def run_pristine_seesaw(fastq_ref_fn, fastq_in1_fn, fastq_in2_fn, fastq_out1_fn, fastq_out2_fn, 
run_record, num_threads, test):
    """run pristine_seesaw.py to ensure paired filtered reads"""

    command = 'pristine_seesaw.py %s %s %s %s %s %d' % (fastq_ref_fn, fastq_in1_fn, fastq_in2_fn, fastq_out1_fn, fastq_out2_fn, num_threads)

    if test:
        print "=== The command ==="
        print command
        return run_record
    
    start = time.time()
    pristine_seesaw.make_pristine_pair(fastq_ref_fn, fastq_in1_fn, fastq_in2_fn, fastq_out1_fn, fastq_out2_fn, num_threads)
    end = time.time()
    run_record.addMessage(program_name=command,
            error_type=LOG_INFO, message='Time taken (mins)',
            value=((end-start)/60.))
    return run_record


def concatenate(pristine_path, contaminated_path, combined_path, run_record, test):
    """concatenates pristine to the end of contaminated producing a new file"""
    def write(infile, outfile):
        n = 0
        for line in infile:
            outfile.write(line)
            n += 1
        return n
    
    total_lines = 0
    concatenated = open(combined_path, 'w')
    
    contaminated = open(contaminated_path, 'r')
    total_lines += write(contaminated, concatenated)
    contaminated.close()
    
    pristine = open(pristine_path, 'r')
    total_lines += write(pristine, concatenated)
    pristine.close()
    
    concatenated.close()
    
    run_record.addMessage(program_name='concatenate',
            error_type=LOG_INFO, message='Total lines', value=total_lines)
    return run_record


def run_bowtie(align_index, fastq_filename, map_filename, run_record, test):
    """run bowtie, output SAM format, add version to run_record"""
    command = 'bowtie --version'
    returncode, stdout, stderr = run_command(command, test)
    for line_ in stdout.splitlines():
        if line_.startswith('bowtie version'):
            run_record.addMessage(program_name=command,
            error_type=LOG_INFO, message='bowtie version: %s' % line_, value=0)

    command = 'bowtie -q --solexa1.3-quals -S --mapq 37 -t -m 1 -p 6 %s %s %s' % \
        (align_index, fastq_filename, map_filename)
    start = time.time()
    returncode, stdout, stderr = run_command(command, test)
    end = time.time()
    
    run_record.addMessage(program_name='bowtie',
            error_type=LOG_INFO, message='Time taken (mins)',
            value=((end-start)/60.))
    
    pipes = {'stderr': stderr, 'stdout': stdout}
    logmsg = {'stderr': LOG_ERROR, 'stdout': LOG_INFO}
    for pipe in pipes:
        for line in pipes[pipe].splitlines():
            print line
            if not line.startswith('#'):
                continue
            line = [element.strip() for element in line[2:].split(':')]
            run_record.addMessage(program_name=command,
                    error_type=logmsg[pipe], message=line[0], value=line[1])
    
    return run_record

def run_bwa_aln(align_index, fastq_filename, out_filename, run_record, num_threads, test):
    """run bwa and add version to run_record"""

    command = 'bwa'
    returncode, stdout, stderr = run_command(command, test)
    for line in stderr.splitlines():
        if line.startswith('Version'):
            run_record.addMessage(program_name=command,
            error_type=LOG_INFO, message=line, value=0) 

    command = 'bwa aln -I -t %s %s %s > %s' % (num_threads, align_index, fastq_filename,
                                        out_filename)
    start = time.time()
    returncode, stdout, stderr = run_command(command, test)
    end = time.time()
    
    run_record.addMessage(program_name=command,
            error_type=LOG_INFO, message='Time taken (mins)',
            value=((end-start)/60.))
    
    pipes = {'stderr': stderr, 'stdout': stdout}
    logmsg = {'stderr': LOG_ERROR, 'stdout': LOG_INFO}
    for pipe in pipes:
        for line in pipes[pipe].splitlines():
            print line
            if not line.startswith('#'):
                continue
            line = [element.strip() for element in line[2:].split(':')]
            run_record.addMessage(program_name=command,
                    error_type=logmsg[pipe], message=line[0], value=line[1])
    
    return run_record

def run_bwa_samse(align_index, sai_filename, fastq_filename,
                    sam_filename, run_record, test):
    """run bwa samse to convert internal coordinates to SAM format output"""
    command = 'bwa samse -n 1 %s %s %s > %s' % (align_index, sai_filename,
            fastq_filename, sam_filename)
    start = time.time()
    returncode, stdout, stderr = run_command(command, test)
    end = time.time()
    
    run_record.addMessage(program_name=command,
            error_type=LOG_INFO, message='Time taken (mins)',
            value=((end-start)/60.))
    
    pipes = {'stderr': stderr, 'stdout': stdout}
    logmsg = {'stderr': LOG_ERROR, 'stdout': LOG_INFO}
    for pipe in pipes:
        for line in pipes[pipe].splitlines():
            print line
            if not line.startswith('#'):
                continue
            line = [element.strip() for element in line[2:].split(':')]
            run_record.addMessage(program_name=command,
                    error_type=logmsg[pipe], message=line[0], value=line[1])
    
    return run_record


def run_bwa_sampe(align_index, in1_sai, in2_sai, in1_fastq, in2_fastq,
            sam_filename, run_record, num_threads, test):
    """runn bwa for paired end reads"""

    command = 'bwa sampe %s %s %s %s %s | pigz -p %s > %s' % (align_index,
                        in1_sai, in2_sai, in1_fastq, in2_fastq, num_threads-1, sam_filename)
    
    start = time.time()
    returncode, stdout, stderr = run_command(command, test)
    end = time.time()
    
    run_record.addMessage(program_name=command,
            error_type=LOG_INFO, message='Time taken (mins)',
            value=((end-start)/60.))
    
    pipes = {'stderr': stderr, 'stdout': stdout}
    logmsg = {'stderr': LOG_ERROR, 'stdout': LOG_INFO}
    for pipe in pipes:
        for line in pipes[pipe].splitlines():
            print line
            if not line.startswith('#'):
                continue
            line = [element.strip() for element in line[2:].split(':')]
            run_record.addMessage(program_name=command,
                    error_type=logmsg[pipe], message=line[0], value=line[1])
    
    return run_record


def convert_sam_to_bam(sam_filename, bam_filename, run_record, test):
    """uses samtools to convert sam to bam"""
    command = 'samtools view -bS %s -o %s' % (sam_filename, bam_filename)
    
    start = time.time()
    returncode, stdout, stderr = run_command(command, test)
    end = time.time()
    
    run_record.addInfo(program_name=command, message='Time taken (mins)',
            value=((end-start)/60.))
    
    pipes = {'stderr': stderr, 'stdout': stdout}
    logmsg = {'stderr': LOG_ERROR, 'stdout': LOG_INFO}
    for pipe in pipes:
        for line in pipes[pipe].splitlines():
            print line
            if not line.startswith('#'):
                continue
            line = [element.strip() for element in line[2:].split(':')]
            run_record.addMessage(program_name=command,
                    error_type=logmsg[pipe], message=line[0], value=line[1])
    
    return run_record

def convert_bam_to_sam(bam_filename, sam_filename, run_record, test):
    """uses samtools to convert sam to bam
        Required for reduce step"""
    command = 'samtools view -h %s -o %s' % (bam_filename, sam_filename)

    start = time.time()
    returncode, stdout, stderr = run_command(command, test)
    end = time.time()

    run_record.addInfo(program_name=command, message='Time taken (mins)',
            value=((end-start)/60.))

    pipes = {'stderr': stderr, 'stdout': stdout}
    logmsg = {'stderr': LOG_ERROR, 'stdout': LOG_INFO}
    for pipe in pipes:
        for line in pipes[pipe].splitlines():
            print line
            if not line.startswith('#'):
                continue
            line = [element.strip() for element in line[2:].split(':')]
            run_record.addMessage(program_name=command,
                    error_type=logmsg[pipe], message=line[0], value=line[1])

    return run_record


def convert_sam_to_sorted_bam(sam_filename, bam_filename, run_record, num_threads, test):
    """uses samtools to convert sam to sorted bam"""
    bam_prefix=bam_filename
    if bam_prefix.endswith('.bam'):
        bam_prefix=bam_prefix[:-4]
    command = 'samtools view -uS %s | samtools sort -m 10000000000 - %s' % (sam_filename, bam_prefix)
    
    start = time.time()
    returncode, stdout, stderr = run_command(command, test)
    end = time.time()
    
    run_record.addInfo(program_name=command, message='Time taken (mins)',
            value=((end-start)/60.))
    
    pipes = {'stderr': stderr, 'stdout': stdout}
    logmsg = {'stderr': LOG_ERROR, 'stdout': LOG_INFO}
    for pipe in pipes:
        for line in pipes[pipe].splitlines():
            print line
            if not line.startswith('#'):
                continue
            line = [element.strip() for element in line[2:].split(':')]
            run_record.addMessage(program_name=command,
                    error_type=logmsg[pipe], message=line[0], value=line[1])
    
    return run_record


def annotate_bam_header(sam_filename, bam_filename, anno_bam_filename, sample_name, run_record, test):
    """write out existing sam header and insert @RG line"""
    sam_prefix=sam_filename
    if sam_prefix.endswith('.sam.gz'):
        sam_prefix=sam_prefix[:-7]
    command = 'samtools view -SH %s > %s.samh' % (sam_filename, sam_prefix)
    file_prefix=sam_prefix.split('/')[-1]
    start = time.time()
    returncode, stdout, stderr = run_command(command, test)
    if sample_name=='': sample_name=file_prefix
    f=open(sam_prefix+'.samh','a')
    f.write("@RG\tPL:ILLUMINA\tSM:%s\tID:%s\n" % (sample_name, file_prefix))
    f.close()
    end = time.time()
    run_record.addInfo(program_name=command, message='Time taken (mins)',
            value=((end-start)/60.))
    
    command = 'samtools reheader %s.samh %s > %s' % (sam_prefix, bam_filename, anno_bam_filename)
    start = time.time()
    returncode, stdout, stderr = run_command(command, test)
    end = time.time()
    run_record.addInfo(program_name=command, message='Time taken (mins)',
            value=((end-start)/60.))
    
    pipes = {'stderr': stderr, 'stdout': stdout}
    logmsg = {'stderr': LOG_ERROR, 'stdout': LOG_INFO}
    for pipe in pipes:
        for line in pipes[pipe].splitlines():
            print line
            if not line.startswith('#'):
                continue
            line = [element.strip() for element in line[2:].split(':')]
            run_record.addMessage(program_name=command,
                    error_type=logmsg[pipe], message=line[0], value=line[1])

    return run_record


def bwa_sampe_to_sorted_bam(align_index, in1_sai, in2_sai, in1_fastq, in2_fastq,
            bam_filename, run_record, sample_name, mem_usage, test):
    """runn bwa for paired end reads"""
    bam_prefix = bam_filename
    if bam_prefix.endswith('.bam'):
        bam_prefix = bam_prefix[:-4]
    dir_split = bam_prefix.split('/')
    sample=sample_name
    if sample=='':
        sample=dir_split[-2]
    rg_line = "@RG\tPL:ILLUMINA\tSM:%s\tID:%s" % (sample, dir_split[-1])
    command = "bwa sampe -r '%s' %s %s %s %s %s | samtools view -buS - | samtools sort -m %d - %s" % (rg_line, align_index, in1_sai, in2_sai, in1_fastq, in2_fastq, mem_usage, bam_prefix)
    # removed -P option which requires 1.25x genome-size RAM to load the entire FM index to reduce disk IO, as this exceeded available memory on nodes
    start = time.time()
    returncode, stdout, stderr = run_command(command, test)
    end = time.time()
    
    run_record.addMessage(program_name=command,
            error_type=LOG_INFO, message='Time taken (mins)',
            value=((end-start)/60.))
    
    pipes = {'stderr': stderr, 'stdout': stdout}
    logmsg = {'stderr': LOG_ERROR, 'stdout': LOG_INFO}
    for pipe in pipes:
        for line in pipes[pipe].splitlines():
            print line
            if not line.startswith('#'):
                continue
            line = [element.strip() for element in line[2:].split(':')]
            run_record.addMessage(program_name=command,
                    error_type=logmsg[pipe], message=line[0], value=line[1])
    
    return run_record

## the method below was to push all things through at the final stage, and probably needs to be reworked into individual steps
def bwa_sampe_to_local_bam_sort_index_move(align_index, in1_sai, in2_sai, in1_fastq, in2_fastq,
            bam_filename, run_record, sample_name, mem_usage, test):
    """runn bwa for paired end reads"""
    bam_prefix = bam_filename
    if bam_prefix.endswith('.bam'):
        bam_prefix=bam_prefix[:-4]
    dir_split = bam_prefix.split('/')
    sample=sample_name
    if sample=='':
        sample=dir_split[-2]
    bam_id=dir_split[-1]
    bam_unsorted_name = "/data/scratch/depressed/" + bam_id + ".unsorted.bam"
    rg_line = "@RG\tPL:ILLUMINA\tSM:%s\tID:%s" % (sample, bam_id)
    command = "bwa sampe -P -r '%s' %s %s %s %s %s | samtools view -b1St %s -o %s -" % \
               (rg_line, align_index, in1_sai, in2_sai, in1_fastq, in2_fastq, align_index+".fai", bam_unsorted_name)
    # N.B. enabling the sampe -P option requires 1.25x genome-size RAM to load the entire FM index to reduce disk IO,
    # and this may exceeded available memory on nodes if more than one job is run at a time
    start = time.time()
    returncode, stdout, stderr = run_command(command, test)
    end = time.time()
    
    run_record.addMessage(program_name=command,
            error_type=LOG_INFO, message='Time taken (mins)',
            value=((end-start)/60.))
    
    pipes = {'stderr': stderr, 'stdout': stdout}
    logmsg = {'stderr': LOG_ERROR, 'stdout': LOG_INFO}
    for pipe in pipes:
        for line in pipes[pipe].splitlines():
            print line
            if not line.startswith('#'):
                continue
            line = [element.strip() for element in line[2:].split(':')]
            run_record.addMessage(program_name=command,
                    error_type=logmsg[pipe], message=line[0], value=line[1])

	local_bam_prefix=bam_unsorted_name[:-13]
    command2 = "samtools sort -m %d %s %s" % (mem_usage, bam_unsorted_name, local_bam_prefix)
    start = time.time()
    returncode, stdout, stderr = run_command(command2, test)
    end = time.time()
    
    run_record.addMessage(program_name=command2,
            error_type=LOG_INFO, message='Time taken (mins)',
            value=((end-start)/60.))
    
    pipes = {'stderr': stderr, 'stdout': stdout}
    logmsg = {'stderr': LOG_ERROR, 'stdout': LOG_INFO}
    for pipe in pipes:
        for line in pipes[pipe].splitlines():
            print line
            if not line.startswith('#'):
                continue
            line = [element.strip() for element in line[2:].split(':')]
            run_record.addMessage(program_name=command2,
                    error_type=logmsg[pipe], message=line[0], value=line[1])
    
    command3 = "samtools index %s" % (local_bam_prefix+'.bam')
    start = time.time()
    returncode, stdout, stderr = run_command(command3, test)
    end = time.time()
    
    run_record.addMessage(program_name=command3,
            error_type=LOG_INFO, message='Time taken (mins)',
            value=((end-start)/60.))
    
    pipes = {'stderr': stderr, 'stdout': stdout}
    logmsg = {'stderr': LOG_ERROR, 'stdout': LOG_INFO}
    for pipe in pipes:
        for line in pipes[pipe].splitlines():
            print line
            if not line.startswith('#'):
                continue
            line = [element.strip() for element in line[2:].split(':')]
            run_record.addMessage(program_name=command3,
                    error_type=logmsg[pipe], message=line[0], value=line[1])

    command4 = "rsync -P %s %s cg:/data/scratch/depressed/%s/%s" % \
        (local_bam_prefix+'.bam', local_bam_prefix+'.bam.bai', sample, bam_id+'.bam')
    start = time.time()
    returncode, stdout, stderr = run_command(command4, test)
    if returncode==0:
        os.remove(local_bam_prefix+'.bam')
        os.remove(local_bam_prefix+'.bam.bai')
    end = time.time()
    
    run_record.addMessage(program_name=command4,
            error_type=LOG_INFO, message='Time taken (mins)',
            value=((end-start)/60.))
    
    pipes = {'stderr': stderr, 'stdout': stdout}
    logmsg = {'stderr': LOG_ERROR, 'stdout': LOG_INFO}
    for pipe in pipes:
        for line in pipes[pipe].splitlines():
            print line
            if not line.startswith('#'):
                continue
            line = [element.strip() for element in line[2:].split(':')]
            run_record.addMessage(program_name=command4,
                    error_type=logmsg[pipe], message=line[0], value=line[1])
    
    return run_record


def index_bam(bam_filename, run_record, test):
    """create .bai index for final .bam file"""
    command = 'samtools index %s' % (bam_filename)
    start = time.time()
    returncode, stdout, stderr = run_command(command, test)
    end = time.time()
    run_record.addInfo(program_name=command, message='Time taken (mins)',
            value=((end-start)/60.))
    
    pipes = {'stderr': stderr, 'stdout': stdout}
    logmsg = {'stderr': LOG_ERROR, 'stdout': LOG_INFO}
    for pipe in pipes:
        for line in pipes[pipe].splitlines():
            print line
            if not line.startswith('#'):
                continue
            line = [element.strip() for element in line[2:].split(':')]
            run_record.addMessage(program_name=command,
                    error_type=logmsg[pipe], message=line[0], value=line[1])

    return run_record


