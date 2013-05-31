"""parses file dumped from R"""
from cogent.parse.table import SeparatorFormatParser, ConvertFields
from cogent import LoadTable
from chippy.util.run_record import RunRecord
from chippy.util.definition import LOG_DEBUG, LOG_INFO, LOG_WARNING, \
    LOG_ERROR, LOG_CRITICAL

__author__ = 'Gavin Huttley, Cameron Jack'
__copyright__ = 'Copyright 2011, Anuj Pahwa, Gavin Huttley, Cameron Jack'
__credits__ = ['Gavin Huttley', 'Cameron Jack']
__license__ = 'GPL'
__maintainer__ = 'Cameron Jack'
__email__ = 'cameron.jack@anu.edu.au'
__status__ = 'pre-release'
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

def validate_exps(probeset_ids, probeset_exps):
    new_ids = []
    new_exps = []

    # takes care of uneven probe/expr info
    # breaks current test for unevenness
    if len(probeset_ids) != len (probeset_exps):
        raise RuntimeError('Different number of probesets to expression values')

    for i in range (len(probeset_ids)):
        if type(probeset_exps[i]) == str:
            continue
        new_ids.append(probeset_ids[i])
        new_exps.append(probeset_exps[i])

    return new_ids, new_exps

def simpleRdumpToTable(path, sep='\t', stable_id_label='', probeset_label='',
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

    Needs to be able to ignore string entries in place of float/double
    """
    rr = RunRecord('simpleRdumpToTable')
    
    # forces reading in as string
    converter = ConvertFields([])
    reader = SeparatorFormatParser(with_title=True, converter=converter,
            sep='\t')
    
    table = LoadTable(path, reader=reader)
    
    # convert '|' delimited probeset and expression fields to tuples of
    # (ints / string) and float respectively
    convert_probeset = convert(to_float=False)
    convert_exp = convert(to_float=True)
    probeset_index = table.Header.index(probeset_label)
    exp_index = table.Header.index(exp_label)
    rows = table.getRawData()

    lost_probesets = 0
    new_rows = []
    for row in rows:
        new_ids = convert_probeset(row[probeset_index])
        new_exps = convert_exp(row[exp_index])
        try:
            new_ids, new_exps = validate_exps(new_ids, new_exps)
        except RuntimeError:
            lost_probesets += 1
            continue
        row[probeset_index] = new_ids
        row[exp_index] = new_exps
        if len(new_ids) == 0:
            continue
        new_rows.append(row)
    
    table = LoadTable(header=table.Header, rows=new_rows)
    if lost_probesets:
        rr.addCritical('probesets not matched to expression values',
                lost_probesets)
    
    if validate:
        assert probeset_label and exp_label and stable_id_label,\
                'Must provide all required column labels to validate'
        stable_ids = table.getDistinctValues(stable_id_label)
        if len(stable_ids) != table.Shape[0]:
            raise RuntimeError('Non unique stable IDs')
        
        rr.addInfo('validation', 'genes were all unique')
        rr.addInfo('validation', 'numbers of probesets matched expression records')

    # look for cases where a probeset maps to multiple genes
    if not allow_probeset_many_gene:
        rr.addInfo('Probesets must map to a single gene', True)
        
        probeset_gene = {}
        gene_problem_probesets = {}
        bad_probesets = set()
        for row in table:
            stable_id = row[stable_id_label]
            for probeset in row[probeset_label]:
                try: # bad probeset, record it and affected genes
                    probeset_gene[probeset].update([stable_id])
                    for gene_id in probeset_gene[probeset]:
                        gene_problem_probesets[gene_id] = \
                                    gene_problem_probesets.get(gene_id, set())
                        gene_problem_probesets[gene_id].update([probeset])
                    
                    bad_probesets.update([probeset])
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
        
        rr.addInfo('No. dropped probesets', len(bad_probesets))
        rr.addInfo('No. genes dropped in probeset filter',
                table.Shape[0]-len(rows))
        table = LoadTable(header=table.Header, rows=rows)

    return table
