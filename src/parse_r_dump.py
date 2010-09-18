"""parses file dumped from R"""
from cogent.parse.table import SeparatorFormatParser, ConvertFields
from cogent import LoadTable

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

def RDumpParser(data):
    """docstring for RDumpParser"""
    if type(data) == str:
        data = open(data)
    
    converter = ConvertFields([(0,long),(1, chrom_num),
                    (2, standardised_strand),
                    (3, long), (4, long), (5, sans_quotes),
                    (6, split_ensembl_ids)])
    parser = SeparatorFormatParser(with_title=True, converter=converter,
            sep='\t')
    reader = parser(data)
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

if __name__ == "__main__":
    infile='../tests/data/rdump_sample.txt'
    table = RDumpToTable(infile)
    print table.getColumns(table.Header[:-1])
