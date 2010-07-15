from cogent import LoadTable
from cogent.parse.table import ConvertFields

# The 4th and the 7th elements of the row of data returned from bowtie are
# integer values and can thus be converted.
row_converter = ConvertFields([(3, int), (6, int)])

def BowtieOutputParser(data):
    """yields a header and row of data from the default bowtie output"""

    header = ['Query Name', 'Strand Direction','Reference Name', 'Offset',
              'Query Seq', 'Quality', 'Other Matches', 'Mismatches']
    yield header

    # If given a filename for the data
    if type(data) == str:
        data = open(data)

    for record in data:
        row = row_converter(record.rstrip('\n').split('\t'))

        # convert the last element to a list of strings
        if row[-1] is '':
            row[-1] = []
        else:
            row[-1] = row[-1].split(',')

        yield row

def BowtieToTable(data):
    """Converts bowtie output to a table"""

    parser = BowtieOutputParser(data)
    header = parser.next()
    rows = [row for row in parser]
    table = LoadTable(header=header, rows=rows)
    return table