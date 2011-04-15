"""parses file dumped from R"""
from cogent.parse.table import SeparatorFormatParser, ConvertFields
from cogent import LoadTable

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"
__version__ = '0.1'

def sans_quotes(val):
    """returns string without quotes"""
    return val.replace("'", "").replace('"', '')

def standardised_strand(val):
    """converts character representation into -1/1 to be consistent with the
    mapping software format"""
    val = sans_quotes(val)
    if val == '+':
        val = 1
    elif val == '-':
        val = -1
    else:
        return None
    
    return val

def split_ensembl_ids(val):
    """returns the ensembl IDs split into elements"""
    val = sans_quotes(val)
    ids = val.split('|')
    return ids

def chrom_num(val):
    """returns just the chromosome ID"""
    val = sans_quotes(val)
    return val.replace('chr', '')

def cast_num(val):
    """int or float"""
    try:
        val = int(val)
    except ValueError:
        val = float(val)
    return val

def get_reader(infile):
    """generates a reader for the infile"""
    for num, line in enumerate(infile):
        if num == 0:
            line = map(sans_quotes, line.split('\t'))
            strand_index = line.index('strand')
            id_index = line.index('ENSEMBL')
        else:
            break
    
    conversions = []
    for index, field in enumerate(line.strip().split('\t')):
        if index == strand_index:
            conversions += [(index, standardised_strand)]
        elif index == id_index:
            conversions += [(index, split_ensembl_ids)]
        else:
            field = eval(field)
            type_ = type(field)
            if type_ == str:
                type_ = sans_quotes
            elif type_ in (int, float):
                type_ = cast_num
        
            conversions += [(index, type_)]
    
    infile.seek(0)
    converter = ConvertFields(conversions)
    reader = SeparatorFormatParser(with_title=True, converter=converter,
            sep='\t')
    return reader(infile)

def RDumpParser(data):
    """docstring for RDumpParser"""
    if type(data) == str:
        data = open(data)
    reader = get_reader(data)
    header = map(sans_quotes, reader.next())
    yield header
    for line in reader:
        yield line

def RDumpToTable(data):
    """returns a cogent table object"""
    parser = RDumpParser(data)
    header = parser.next()
    rows = [l for l in parser]
    return LoadTable(header=header, rows=rows, space=2)

def convert(to_float=False, strict=False):
    """converts a | separated string into a tuple of floats or ints"""
    type_ = [int, float][to_float]
    def call(value):
        result = value.strip().split('|')
        try:
            result = tuple(map(type_, result))
        except ValueError:
            # can occur as older Affy arrays have probeset IDs that mix
            # int/str
            if strict:
                raise ValueError
            pass
        
        return result
    
    return call

def SimpleRdumpToTable(path, sep='\t', stable_id_label='', probeset_label='', exp_label='', validate=True):
    """returns a cogent table object
    
    Handles case where probset id's and expressions scores are separated by
    the pipe -- | -- character. The probset and expression scores are then
    converted to tuples of ints or floats respectively.
    
    Arguments:
        - probeset_label: name of column containing probesets
        - exp_label: name of column containing expression scores
        - stable_id_label: name of column containing Ensembl stable IDs
        - validate: checks that -- stable IDs are unique in the file, that
          for each row the number of probesets equals the number of expression
          scores. Raises a RuntimeException if failure occurs for any
          of these checks.
    """
    
    converter = ConvertFields([(1, convert(to_float=False)),
                               (2, convert(to_float=True))])
    reader = SeparatorFormatParser(converter=converter, sep='\t')
    table = LoadTable(path, reader=reader)
    if validate:
        assert probeset_label and exp_label and stable_id_label,\
            'Must provide all required column labels to validate'
        stable_ids = table.getDistinctValues(stable_id_label)
        if len(stable_ids) != table.Shape[0]:
            raise RuntimeError('Non unique stable IDs')
        
        for row in table:
            if len(row[probeset_label]) != len(row[exp_label]):
                raise RuntimeError(
                        'Mismatched number of probesets and exp scores')
    
    return table

