import re
from cogent.parse.table import ConvertFields, SeparatorFormatParser

__author__ = "Anuj Pahwa, Gavin Huttley, Cameron Jack"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley, Cameron Jack"
__credits__ = ["Anuj Pahwa", "Gavin Huttley, Cameron Jack"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"
__version__ = '0.1'

_not_digits = re.compile(r'\D+')

def _strict_cigar_span(val):
    return sum([int(v) for v in _not_digits.split(val) if v])

def _cigar_span(val):
    try:
        r = int(val[:-1])
    except ValueError:
        r = 0
    return r

def get_strand(val):
    """returns 1/-1 for strand from bitwise operation"""
    v = int(val)
    strand = [-1,1][v & 16 == 0]
    return strand

def zero_based(val):
    """returns a zero-based integer"""
    return int(val) - 1

strict_converter = ConvertFields([(1, int),(3,int),(4,int),
                                    (5, _strict_cigar_span)])

converter = ConvertFields([(1, get_strand), (3, zero_based),(4,int), (5, _cigar_span)])

# SAM fields: QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL, OPT
complete_converter = ConvertFields([(0, str), (1, get_strand), (2, str), (3, zero_based),
                                    (4, int), (5, _cigar_span), (6, str), (7, int),
                                    (8, int), (9, str), (10, str), (11, str)])

def MinimalSamParser(data, converter=converter):
    """returns records from a sam file
    
    NOTE: the default converter turns the 1-based numbering of POS into
    0-based numbering"""
    # If given a filename for the data
    if type(data) == str:
        data = open(data)
    
    # get the lengths dict
    lengths = {}
    for row in data:
        if not row.startswith('@'):
            yield lengths
            break
        elif not row.startswith('@SQ'):
            continue
        line = row.split()[1:]
        name = line[0].split(':')[1]
        length = int(line[1].split(':')[1])
        lengths[name] = length
    
    parser = SeparatorFormatParser(converter=converter, with_header=False,
                                   sep="\t")
    
    for row in parser(data):
        yield row

def CompleteSamParser(data, converter=complete_converter):
    """returns records from a sam file

    NOTE: the default converter turns the 1-based numbering of POS into
    0-based numbering"""
    # If given a filename for the data
    if type(data) == str:
        data = open(data)

    # get the lengths dict
    lengths = {}
    header_lines = 0

    for row in data:
        header_lines += 1
        if not row.startswith('@'):
            yield lengths
            break
        elif not row.startswith('@SQ'):
            continue
        line = row.split()[1:]
        name = line[0].split(':')[1]
        length = int(line[1].split(':')[1])
        lengths[name] = length

    data.seek(0)
    for i, line in enumerate(data):
        if i == header_lines - 2:
            break

    parser = SeparatorFormatParser(converter=complete_converter, with_header=False,
                                   sep="\t")

    for row in parser(data):
        yield row
