import time
from cogent import LoadTable

from parse_bowtie import BowtieOutputParser
from parse_sam import MinimalSamParser
from region_count import RegionCounts

def format_long(n):
    """format a long integer"""
    n = str(int(n))
    num, rem = divmod(len(n),3)
    r = []
    if rem:
        r = [n[:rem]]

    for i in range(rem, len(n), 3):
        r += [n[i: i+3]]

    return ','.join(r)

# following needs to be moved to client code
def bufcount(filename):
    # from http://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python
    f = open(filename)
    lines = 0
    buf_size = 1024 * 1024
    read_f = f.read # loop optimization

    buf = read_f(buf_size)
    while buf:
        lines += buf.count('\n')
        buf = read_f(buf_size)

    return lines

def get_file_length(filename):
    """output number of lines in read file"""

    start = time.time()
    print 'Reading %s' % filename
    num_lines = bufcount(filename)
    end = time.time()
    print format_long(num_lines) + ' lines read in ' + \
            ('%fs' % ((end-start)/60.))

def get_read_counts_sam(infile_name, chrom_name):
    """Generate the read counts using the SAM format infile"""

    parser = MinimalSamParser(infile_name)
    n = 0
    lengths = parser.next()
    counter = RegionCounts(lengths[chrom_name])
    print 'Parsing: %s' % infile_name
    for record in parser:
        if n % 100000 == 0:
            print format_long(n)
        n += 1
        chrom = record[2]
        if chrom != chrom_name:
            continue

        start = record[3]
        span = record[5]

        if span:
            counter.addRead(start, start+span)

    print '\nTotal %s counts: %s' % (chrom_name,
                            counter.counts.sum())
    print 'Num reads: %d' % n
    return counter

def get_read_counts_bowtie(infile_name, chrom_name, chrom_size):
    """Generate the read counts using the bowtie format infile"""

    parser = BowtieOutputParser(infile_name)
    header = parser.next()
    counter = RegionCounts(chrom_size, one_based=False)
    n = 0
    for record in parser:
        n += 1
        chrom = record[2]
        if chrom != chrom_name:
            continue

        start = record[3]
        span = len(record[4])

        if span:
            counter.addRead(start, start+span)

    print '\nTotal %s counts: %s' % (chrom_name,
                            counter.counts.sum())
    print 'Num reads: %d' % n
    return counter

if __name__ == "__main__":

    # Example of calling the sam method
    # get_read_counts_sam('../data/s_8_pristine.map', 'chr1')

    # manually supplying the length of the chromosome for testing purposes.
    get_file_length('../data/bowtie_output.map')
    get_read_counts_bowtie('../data/bowtie_output.map', 'chr1', 197195432)
