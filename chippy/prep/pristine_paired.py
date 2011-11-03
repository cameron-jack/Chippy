#!/usr/bin/env python26
import sys
import time
from subprocess import Popen, PIPE

from cogent.parse.fastq import MinimalFastqParser

def run_command(command):
    """executes a command"""
    r = Popen(command, shell=True, stdout=PIPE, stderr=PIPE,
                bufsize=1, close_fds=True)
    return r

def get_out_pipe(command):
    """opens a stdin pipe"""
    r = Popen(command, shell=True, stdin=PIPE, bufsize=1, close_fds=True)
    return r

def load_labels(input_fn, num_procs):
    """uses grep and cut to build dict"""
    command = "pigz -dcp %(dthreads)s %(infile)s " % dict(dthreads=num_procs,
                                                    infile=input_fn)
    command += "| parallel --pipe 'grep @ | cut -d: -f3,4,5 | cut -d# -f1'"
    
    p = run_command(command)
    start = time.time()
    # critical to use a generator NOT list comprehension here
    labels = dict((l.strip(), 0) for l in p.stdout)
    end = time.time()
    p.stdout.close()
    p.stderr.close()
    return labels

#def get_start_label_slice(label):
#    """returns unique component of fastq sequence label"""
#    split = label.split(':')
#    return len(split[0]) + len(split[1]) + 2

def get_unique_portion(label):
    """returns unique component of fastq sequence label (all but 1st 2 colon-delim ints)"""
# to include handling for CASAVA 1.8+
#    return label.split(':',2)[2].rsplit('#',1)[0].rsplit(' ',1)[0] 
    return label.split(':',2)[2].rsplit('#',1)[0]

def write_fastq(out, label, seq, quals):
    """write fastq formatted seq record"""
    out.write("@%s\n%s\n+\n%s\n" % (label, seq, quals))

def parse_fastq(input_fn, output_fn, labels, num_procs, set_val):
    """parses a fastq file only writing reads in labels if set_val is True
    or only writing reads in labels with value=1 if set_val is False"""
    # rip the file using pigz
    
    dprocs = 1
    if num_procs > 7:
        dprocs=2
    cprocs = num_procs - dprocs
    if cprocs < 1:
        cprocs=1    
    decomp_command = "pigz -dcp %(dthreads)s %(infile)s " % dict(
        dthreads=dprocs, infile=input_fn)
    decompress = run_command(decomp_command)
    parser = MinimalFastqParser(decompress.stdout, strict=False)
    comp_command = "pigz -cp %(cthreads)s > %(outfile)s " % dict(
        cthreads=cprocs, outfile=output_fn)
    compress = get_out_pipe(comp_command)
#    start_slice = 1
    num_reads_tot = 0
    num_reads_write = 0
    if set_val:
        for name, seq, qual in parser:
            num_reads_tot += 1
            uniq = get_unique_portion(name)
#            if name[start_slice-1] != ':': # allow possible diff length starts
#                start_slice = get_start_label_slice(name)
#        
            if uniq not in labels:
                continue
            
            # write out data here
            labels[uniq] = 1
            num_reads_write += 1
            write_fastq(compress.stdin, name, seq, qual)
    else:
        for name, seq, qual in parser:
            num_reads_tot += 1
            uniq = get_unique_portion(name)
#            if name[start_slice-1] != ':': # allow possible diff length starts
#                start_slice = get_start_label_slice(name)
            
            if uniq not in labels:
                continue
            if labels[uniq] != 1:
                continue
            
            num_reads_write += 1
            write_fastq(compress.stdin, name, seq, qual)
        
        # we write out result
    decompress.stdout.close()
    decompress.stderr.close()
    compress.stdin.close()
    print 'Total num reads = %d' % (num_reads_tot)
    print 'Write num reads = %d' % (num_reads_write)

def make_pristine_pair(input_fn1, input_fn2, output_fn1, output_fn2, num_procs=8):
    """converts paired input fastqs to pristine paired output"""
    labels = load_labels(input_fn2, num_procs)
    parse_fastq(input_fn1, output_fn1, labels, num_procs, True)
    parse_fastq(input_fn2, output_fn2, labels, num_procs, False)

def usage():
    print 'Usage: pristine_paired.py <input1.fq.gz> <input2.fq.gz> <output1.fq.gz> <output2.fq.gz> <ncpus>'
    sys.exit(-1)

def main(*args):
    try:
        make_pristine_pair(args[1],args[2],args[3],args[4],int(args[5]))
    except:
        usage()
    else:
        return 0 # exit errorlessly
 
if __name__ == '__main__':
    sys.exit(main(*sys.argv))
## 
