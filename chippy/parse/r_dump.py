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

def convert(to_float=False, strict=False):
    """converts a | separated string into a tuple of floats or ints"""
    type_ = [int, float][to_float]
    def call(value):
        result = value.strip().split('|')
        try:
            result = map(type_, result)
        except ValueError:
            # can occur as older Affy arrays have probeset IDs that mix
            # int/str
            if strict:
                raise ValueError
            pass
        
        return tuple(result)
    
    return call

def remove_probesets(row, probesets, probeset_index, exp_index):
    """returns a new row without the problem probesets or their expression"""
    row = list(row)
    current_probesets = list(row[probeset_index])
    current_exp = list(row[exp_index])
    del_indices = map(current_probesets.index, probesets)
    new_probesets = []
    new_exp = []
    for i in range(len(current_probesets)):
        if i in del_indices:
            continue
        
        new_probesets.append(current_probesets[i])
        new_exp.append(current_exp[i])
    
    if len(new_probesets) == 0:
        return None
    
    row[probeset_index] = tuple(new_probesets)
    row[exp_index] = tuple(new_exp)
    return row

def SimpleRdumpToTable(path, sep='\t', stable_id_label='', probeset_label='',
        exp_label='', allow_probeset_many_gene=False, validate=True):
    """returns a cogent table object
    
    Handles case where probset id's and expressions scores are separated by
    the pipe -- | -- character. The probset and expression scores are then
    converted to tuples of ints or floats respectively.
    
    Arguments:
        - probeset_label: name of column containing probesets
        - exp_label: name of column containing expression scores
        - stable_id_label: name of column containing Ensembl stable IDs
        - allow_probeset_many_gene: whether one probeset can map to multiple
          genes
        - validate: checks that -- stable IDs are unique in the file, that
          for each row the number of probesets equals the number of expression
          scores. Raises a RuntimeException if failure occurs for any
          of these checks.
    """
    # forces reading in as string
    converter = ConvertFields([])
    reader = SeparatorFormatParser(with_title=True, converter=converter,
            sep='\t')
    
    table = LoadTable(path, reader=reader)
    
    # convert '|' delimited probeset and expression fields to tuples of
    # (ints / string) and int respectively
    convert_probeset = convert(to_float=False)
    convert_exp = convert(to_float=True)
    probeset_index = table.Header.index(probeset_label)
    exp_index = table.Header.index(exp_label)
    rows = table.getRawData()
    for row in rows:
        row[probeset_index] = convert_probeset(row[probeset_index])
        row[exp_index] = convert_exp(row[exp_index])
    
    table = LoadTable(header=table.Header, rows=rows)
    
    if validate:
        assert probeset_label and exp_label and stable_id_label,\
            'Must provide all required column labels to validate'
        stable_ids = table.getDistinctValues(stable_id_label)
        if len(stable_ids) != table.Shape[0]:
            raise RuntimeError('Non unique stable IDs')
        
        for row in table:
            if len(row[probeset_label]) != len(row[exp_label]):
                gene_id = row[stable_id_label]
                raise RuntimeError('Mismatched number of probesets and '\
                        'exp scores for %s' % gene_id)
            
        
    
    # look for cases where a probeset maps to multiple genes
    if not allow_probeset_many_gene:
        probeset_gene = {}
        gene_problem_probesets = {}
        for row in table:
            stable_id = row[stable_id_label]
            for probeset in row[probeset_label]:
                try: # bad probeset, record it and affected genes
                    probeset_gene[probeset].update([stable_id])
                    for gene_id in probeset_gene[probeset]:
                        gene_problem_probesets[gene_id] = \
                                    gene_problem_probesets.get(gene_id, set())
                        gene_problem_probesets[gene_id].update([probeset])
                
                except KeyError:
                    probeset_gene[probeset] = set([stable_id])
        
        # what indices for the essential columns
        id_index = table.Header.index(stable_id_label)
        probeset_index = table.Header.index(probeset_label)
        exp_index = table.Header.index(exp_label)
        rows = []
        for row in table.getRawData():
            stable_id = row[id_index]
            if stable_id in gene_problem_probesets:
                row = remove_probesets(row, gene_problem_probesets[stable_id],
                            probeset_index, exp_index)
            if row is not None:
                rows.append(row)
        
        table = LoadTable(header=table.Header, rows=rows)
    
    return table

