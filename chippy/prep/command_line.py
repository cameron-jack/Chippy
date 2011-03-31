#!/usr/bin/env python
import os, sys, re, time
import subprocess
from optparse import OptionParser

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley"
__credits__ = ["Gavin Huttley"]
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
    return_value = r.returncode
    return returncode, stdout, stderr

def run_blat(blat_adapters, query_file, psl_out, run_record, test):
    """run blat"""
    command = "blat %s %s %s" % (blat_adapters, query_file, psl_out)
    if test:
        print "=== The command ==="
        print command
        return
    start = time.time()
    returncode, stdout, stderr = run_command(command, test)
    end = time.time()
    run_record.addMessage(program_name='blat',
            error_type='stdout', message='Time taken (mins)',
            value= ((end-start)/60.))
    if stdout:
        print
        print ''.join(stdout)
    
    if stderr:
        print
        print ''.join(stderr)
    return run_record

def run_bowtie(bowtie_index, save_dir, fastq_filename, map_filename, run_record, test):
    """run bowtie"""
    command = 'bowtie -q --solexa1.3-quals -t -m 1 -p 8 %s %s %s' % \
        (bowtie_index, fastq_filename, map_filename)
    start = time.time()
    returncode, stdout, stderr = run_command(command, test)
    end = time.time()
    
    run_record.addMessage(program_name='bowtie',
            error_type='stdout', message='Time taken (mins)',
            value= ((end-start)/60.))
    
    pipes = {'stderr': stderr, 'stdout': stdout}
    for pipe in pipes:
        for line in pipes[pipe].splitlines():
            print line
            if not line.startswith('#'):
                continue
            line = [element.strip() for element in line[2:].split(':')]
            run_record.addMessage(program_name='bowtie',
                        error_type=pipe, message=line[0], value=line[1])
    
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
                error_type='stdout', message='Total lines', value=total_lines)
    return run_record

