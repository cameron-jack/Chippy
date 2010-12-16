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
    get_read_counts_bowtie('../data/bowtie_output.map', 'chr1', 197195432)
