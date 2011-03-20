#!/usr/bin/env python
import os, sys, re
import subprocess
from optparse import OptionParser

def run_command(command, test):
    """executes a command"""
    pipe = subprocess.PIPE
    r = subprocess.Popen(command, shell=True, stdout=pipe, stderr=pipe, bufsize=-1)
    stdout = []
    if r.stdout:
        stdout = r.stdout.readlines()
    stderr = []
    if r.stderr:
        stderr = r.stderr.readlines()
    
    return stdout, stderr

def run_blat(blat_adapters, query_file, psl_out, run_record, test):
    """run blat"""
    command = "blat %s %s %s" % (blat_adapters, query_file, psl_out)
    if test:
        print "=== The command ==="
        print command
        return
    stdout, stderr = run_command(command, test)
    if stdout:
        print
        print ''.join(stdout)
    
    if stderr:
        print
        print ''.join(stderr)
    return run_record

def run_bowtie(bowtie_index, save_dir, fastq_filename, map_filename, run_record, test):
    """run bowtie"""
    fastq_filename = os.path.join(save_dir, fastq_filename)
    map_filename = os.path.join(save_dir, map_filename)
    command = 'bowtie -q --solexa1.3-quals -t -m 1 -p 8 %s %s %s' % \
        (bowtie_index, fastq_filename, map_filename)
    stdout, stderr = run_command(command, test)
    pipes = {'stderr': stderr, 'stdout': stdout}
    for pipe in pipes:
        for line in pipes[pipe]:
            print line
            if not line.startswith('#'):
                continue
            line = [element.strip() for element in line[2:].split(':')]
            run_record.addMessage(program_name='bowtie',
                        error_type=pipe, message=line[0], value=line[1])
    
    return run_record


class RunRecord(object):
    """object for recording program messages"""
    def __init__(self):
        super(RunRecord, self).__init__()
        self.records = []
        
    def addMessage(self, program_name, error_type, message, value):
        """add a message about an execution"""
        self.records.append([program_name, error_type, message, value])
    
    def getMessageTable(self):
        """docstring for display"""
        header = ['program_name', 'error_type', 'message', 'value']
        table = LoadTable(header=header, rows=self.records)
        return table
    
    def display(self):
        table = self.getMessageTable()
        print table
    
