import re
from cogent.parse.table import ConvertFields, SeparatorFormatParser

_not_digits = re.compile(r'\D+')

def _strict_cigar_span(val):
    return sum([int(v) for v in _not_digits.split(val) if v])

def _cigar_span(val):
    try:
        r = int(val[:-1])
    except ValueError:
        r = 0
    return r

strict_converter = ConvertFields([(1, int),(3,int),(4,int),
                                    (5, _strict_cigar_span)])

converter = ConvertFields([(3,int),(4,int), (5, _cigar_span)])

def MinimalSamParser(data):
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
