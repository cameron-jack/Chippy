#!/usr/bin/env python
import os, sys, re, time
import subprocess

from chippy.util.definition import LOG_DEBUG, LOG_INFO, LOG_WARNING, \
    LOG_ERROR, LOG_CRITICAL

__author__ = "Gavin Huttley, Cameron Jack"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley, Cameron Jack"
__credits__ = ["Gavin Huttley", "Cameron Jack"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"
__version__ = '0.1'


def run_command(command, test):
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

def run_fastx_clipper(blat_adapters, fastq_in_fn, fastq_out_fn, run_record, test):
    """run fastx_clipper to remove adapter sequences"""

    file_in = open (blat_adapters)
    # ignore header
    line = file_in.readline()
    # get first sequence
    line = file_in.readline()
    line = line.rstrip('\n')
    line = line.rstrip('\r')
    file_in.close()

    command = "fastx_clipper -a %s -i %s -o %s" % (line, fastq_in_fn, fastq_out_fn)
    if test:
        print "=== The command ==="
        print command
        return
    start = time.time()
    returncode, stdout, stderr = run_command(command, test)
    end = time.time()
    run_record.addMessage(program_name='fastx_clipper',
            error_type=LOG_INFO, message='Time taken (mins)',
            value=((end-start)/60.))
    if stdout:
        print
        print ''.join(stdout)

    if stderr:
        print
        print ''.join(stderr)
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

def run_bwa_aln(align_index, fastq_filename, out_filename, run_record, test):
    """run bwa and add version to run_record"""

    command = 'bwa'
    returncode, stdout, stderr = run_command(command, test)
    for line in stderr.splitlines():
        if line.startswith('Version'):
            run_record.addMessage(program_name=command,
            error_type=LOG_INFO, message=line, value=0) 

    command = 'bwa aln -t 6 %s %s > %s' % (align_index, fastq_filename,
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

