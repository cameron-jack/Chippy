import time
from cogent import LoadTable

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

def get_read_counts_sam(infile_name, chrom_name):
    """Generate the read counts using the SAM format infile"""

    from parse_sam import MinimalSamParser

    start = time.time()
    num_lines = bufcount(infile_name)
    end = time.time()
    print format_long(num_lines)
    print (end-start)/60.
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

        #flag = [0, 1][int(record[1]) & 16 == 16] # 1 if reverse
        start = record[3]-1 #1 based indexing in sam file
        span = record[5]

        if span:
            counter.addRead(start, start+span)

    counter.save(infile_name.replace('.map', '-%s' % chrom_name))
    print '\n\nTotal %s counts: %s' % (chrom_name,
                            counter.counts.sum())
    print 'Num reads: %d' % n

def get_read_counts_bowtie(infile_name, chrom_name, chrom_size):
    """Generate the read counts using the bowtie format infile"""

    from parse_bowtie import BowtieOutputParser

    start = time.time()
    num_lines = bufcount(infile_name)
    end = time.time()
    print format_long(num_lines) + ' lines read in ' + ('%fs' % ((end-start)/60.))
    parser = BowtieOutputParser(infile_name)
    header = parser.next()
    counter = RegionCounts(chrom_size)
    n = 0
    for record in parser:
        n += 1
        if n % 100000 == 0:
            print format_long(n)

        chrom = record[2]
        if chrom != chrom_name:
            continue

        #flag = [0, 1][record[1] == '-'] # 1 if reverse
        start = record[3]
        span = len(record[4])

        if span:
            counter.addRead(start, start+span)

    counter.save(infile_name.replace('.map', '-%s' % chrom_name))
    print '\n\nTotal %s counts: %s' % (chrom_name,
                            counter.counts.sum())
    print 'Num reads: %d' % n



if __name__ == "__main__":
    # get_read_counts('../data/s_7_contaminated.map', 'chr1')
    # get_read_counts('../data/s_8_pristine.map', 'chr1')
    # manually supplying the length of the chromosome for testing purposes.
    get_read_counts_bowtie('../data/bowtie_output.map', 'chr1', 197195432)
