""" parses expr/diff score tab-delimited files """
from cogent.util.table import Table
from chippy.util.run_record import RunRecord
from gzip import GzipFile
from chippy.express.definition import EXPR_HEADER, DIFF_HEADER

__author__ = 'Cameron Jack'
__copyright__ = 'Copyright 2011-2013, Gavin Huttley, Cameron Jack, Anuj Pahwa'
__credits__ = ['Gavin Huttley', 'Cameron Jack']
__license__ = 'GPL'
__maintainer__ = 'Cameron Jack'
__email__ = 'cameron.jack@anu.edu.au'
__status__ = 'pre-release'
__version__ = '0.2'

def _check_expr_headers(header_row, stable_id_label='', probeset_label='',
        exp_label=''):
    """
        Check the header labels match for standard expression. Probeset is
        optional and results in probes_present being False
    """
    rr = RunRecord('_check_expr_headers')

    try:
        gene_col = header_row.index(stable_id_label)
    except ValueError:
        rr.dieOnCritical('Stable ID column header not found in', header_row)

    try:
        exp_col = header_row.index(exp_label)
    except ValueError:
        rr.dieOnCritical('Expression score column header not found in',
                header_row)

    try:
        probe_col = header_row.index(probeset_label)
        probes_present = True
    except ValueError:
        rr.addWarning('Probeset column header not found in', header_row)
        probe_col = -1
        probes_present = False

    return gene_col, probe_col, exp_col, probes_present

def _check_diff_headers(header_row, sig_label='', pval_label=''):
    """ check the additional significance header label and
        p_value label. Has be able to deal with there having been
        no probesets in column 2.
    """
    rr = RunRecord('_check_diff_headers')
    try:
        sig_col = header_row.index(sig_label)
    except ValueError:
        rr.dieOnCritical('Significance column header not found in', header_row)

    try:
        pval_col = header_row.index(pval_label)
    except ValueError:
        rr.addCritical('Expected to see P-value column header', pval_label)
        rr.dieOnCritical('P-value column header not found in', header_row)

    return sig_col, pval_col

def _validate_probes_scores(genes, probes, exp, sig=None, pval=None):
    """
        The number of probes must match the number of scores for each gene.
        This is shared by expr and diff.
    """
    rr = RunRecord('_validate_probes_scores')
    genes_removed = 0
    valid_genes = []
    valid_probes = []
    valid_exp = []
    valid_sigs = []
    valid_pvals = []
    if sig is not None and pval is not None:
        for g, p, e, s, v in zip(genes, probes, exp, sig, pval):
            if len(p) == len(e):
                valid_genes.append(g)
                valid_probes.append(p)
                valid_exp.append(e)
                valid_sigs.append(s)
                valid_pvals.append(v)
            else:
                genes_removed += 1
    else:
        for g, p, e in zip(genes, probes, exp):
            if len(p) == len(e):
                valid_genes.append(g)
                valid_probes.append(p)
                valid_exp.append(e)
            else:
                genes_removed += 1

    rr.addInfo('Number of genes removed in validation step',
            genes_removed)
    if len(valid_sigs) > 0:
        return valid_genes, valid_probes, valid_exp, valid_sigs, valid_pvals

    return valid_genes, valid_probes, valid_exp

def _remove_multimapped_probesets(genes, probes, exp, sig=None, pval=None):
    """
        If a probe maps to more than one gene then remove it and any associated
        score. First, using a dictionary of probe to gene mappings, build sets
        of shared probes and affected genes. Second, go through each gene and
        if affected, go through the probes and keep non-duplicated probes and
        their associated scores. If the gene has any probes/scores remaining
        after this, add it and the probes & scores to the verified lists.

        Shared by expr and diff, so need to be able to deal with sig and
        pval if also present.
    """
    rr = RunRecord('_remove_multimapped_probesets')

    probe_to_gene = {}
    shared_probes = set()
    affected_genes = set()
    for g, p in zip(genes, probes):
        for probe in p:
            if probe in probe_to_gene.keys():
                shared_probes.add(probe)
                affected_genes.add(g)
                affected_genes.add(probe_to_gene[probe])
            else:
                probe_to_gene[probe] = g

    v_genes = []; v_probes = []; v_exp = []; v_sig = []; v_pval = []
    if sig is not None and pval is not None:
        for g, p, e, s, v in zip(genes, probes, exp, sig, pval):
            if g not in affected_genes:
                v_genes.append(g)
                v_probes.append(p)
                v_exp.append(e)
                v_sig.append(s)
                v_pval.append(v)
            elif g in affected_genes:
                nondupe_probes = []
                nondupe_exp = []
                for probe, exp in zip(p, e):
                    if probe not in shared_probes:
                        nondupe_probes.append(probe)
                        nondupe_exp.append(exp)
                if len(nondupe_probes) > 0:
                    v_genes.append(g)
                    v_probes.append(nondupe_probes)
                    v_exp.append(nondupe_exp)
                    v_sig.append(s)
                    v_pval.append(v)
    else:
        for g, p, e in zip(genes, probes, exp):
            if g not in affected_genes:
                v_genes.append(g)
                v_probes.append(p)
                v_exp.append(e)
            elif g in affected_genes:
                nondupe_probes = []
                nondupe_exp = []
                for probe, exp in zip(p, e):
                    if probe not in shared_probes:
                        nondupe_probes.append(probe)
                        nondupe_exp.append(exp)
                if len(nondupe_probes) > 0:
                    v_genes.append(g)
                    v_probes.append(nondupe_probes)
                    v_exp.append(nondupe_exp)
            # or just pass it by if no probes/scores left

    rr.addInfo('Number of genes removed in validation step',
            len(genes)-len(v_genes))
    rr.addInfo('Number of probes/scores removed',
            len(probes)-len(v_probes))

    if len(v_sig) > 0:
        return v_genes, v_probes, v_exp, v_sig, v_pval
    return v_genes, v_probes, v_exp

def _read_data_file(data_path, sep='\t', stable_id_label='ENSEMBL',
        probeset_label='probeset', exp_label='exp', sig_label='sig',
        pval_label='p_val', is_diff=False):
    """
        Get the data out of the file - shared by both:
        gene_expr_to_table() and gene_expr_diff_to_table()
    """
    rr = RunRecord('_read_data_file')

    rows = []
    if '.gz' in data_path.lower():
        with GzipFile(data_path, 'r') as data_file:
            for data in data_file:
                rows.append(str(data).strip().split(sep))
    else:
        with open(data_path) as data_file:
            for data in data_file:
                rows.append(str(data).strip().split(sep))

    # check that headers are valid and whether probesets are present
    gene_col, probe_col, exp_col, probes_present = _check_expr_headers(rows[0],
            stable_id_label=stable_id_label, probeset_label=probeset_label,
            exp_label=exp_label)
    if not probes_present:
        rr.addInfo('No probeset header found. Reading as', 'RNA-seq')
    else:
        rr.addInfo('Probesets found. Reading as', 'Micro-array')

    if is_diff:
        sig_col, pval_col = _check_diff_headers(rows[0], sig_label=sig_label,
                pval_label=pval_label)

    genes = []; probes = []; exp = []; sigs = []; pvals = []

    # get data from each row
    for i, row in enumerate(rows):
        if i==0:
            continue # skip header line

        genes.append(str(row[gene_col]))

        if probes_present:
            probes.append(list(row[probe_col].split('|')))
        else:
            probes.append('P'+str(i)) # give RNAseq data a fake probe id

        exp_strs = row[exp_col].split('|')
        if is_diff:
            sigs.append(str(row[sig_col]))
            pvals.append(float(row[pval_col]))

        # Nuke any scores and probes marked with 'NA' by R
        while 'NA' in exp_strs:
            exp_strs.remove('NA')
        while 'NA' in probes:
            probes.remove('NA')

        if len(exp_strs) == 0:
            rr.dieOnCritical('No expression scores remaining', row)
        if len(probes) == 0:
            rr.dieOnCritical('No probes remaining', row)

        try:
            expression = map(float, exp_strs)
            exp.append(expression)
        except ValueError:
            rr.addCritical('Expected expression score float on line', i)
            rr.dieOnCritical('Line has incorrect format', row)

    if is_diff:
        return genes, probes, exp, sigs, pvals, probes_present
    else:
        return genes, probes, exp, probes_present

def gene_expr_to_table(data_path, sep='\t', stable_id_label='',
        probeset_label='', exp_label='', allow_probeset_many_gene=False,
        validate=True):
    """
        Returns a cogent table object

        Deals with a simple tab-delimited representation of gene expression
        data which may have come from either micro-array or mRNA-seq
        experiments.

        Data from micro-arrays will have probeset information for each
        gene and a score to match each probe.

        RNA-seq data will not have probes and simply a single score for each
        gene. In this case we will create a fake probe for each gene of the
        form 'P' + a unique integer.

        Probset id's and expressions scores are separated by the pipe
        -- | -- character. The probset and expression scores are then
        converted to tuples of ints or floats respectively.

        Arguments:
            - probeset_label: name of column containing probesets
            - exp_label: name of column containing expression scores
            - stable_id_label: name of column containing Ensembl stable IDs
            - allow_probeset_many_gene: whether one probeset can map to
                multiple genes. If not we remove probes and scores that multi-
                map.
            - validate: checks that -- stable IDs are unique in the file,
                that for each row the number of probesets equals the
                number of expression scores. Removes the gene entry.
    """

    rr = RunRecord('geneExprDataToTable')

    rr.addInfo('Reading expression data', data_path)
    genes, probes, exp, probes_present = _read_data_file(data_path, sep=sep,
            stable_id_label=stable_id_label, probeset_label=probeset_label,
            exp_label=exp_label)

    if probes_present:
        if validate:
            # if probes and scores are mismatched, nuke the gene
            genes, probes, exp = \
                    _validate_probes_scores(genes, probes, exp)

        if not allow_probeset_many_gene:
            # each probe should map to only one gene
            genes, probes, exp = \
                    _remove_multimapped_probesets(genes, probes, exp)

    rows = [[g,p,e] for g,p,e in zip(genes, probes, exp)]
    return Table(header=EXPR_HEADER, rows=rows)

def gene_expr_diff_to_table(data_path, sep='\t', stable_id_label='',
        probeset_label='', exp_label='', sig_label='', pval_label='',
        allow_probeset_many_gene=False, validate=True):
    """
        As per gene_expr_to_table() but with the addition of sig_label and
        pval_label columns.
    """
    rr = RunRecord('gene_expr_diff_to_table')

    rr.addInfo('Reading expression diff file', data_path)
    genes, probes, exp, sig, pval, probes_present = _read_data_file(\
            data_path, sep=sep, stable_id_label=stable_id_label,
            probeset_label=probeset_label, exp_label=exp_label,
            sig_label=sig_label, pval_label=pval_label, is_diff=True)

    if probes_present:
        if validate:
            # if probes and exp are mismatched, nuke the gene
            genes, probes, exp, sig, pval =\
                    _validate_probes_scores(genes, probes, exp, sig, pval)

        if not allow_probeset_many_gene:
            # each probe should map to only one gene
            genes, probes, exp, sig, pval =\
                    _remove_multimapped_probesets(genes, probes, exp,
                    sig, pval)

    header = DIFF_HEADER
    rows = [[g, p, e, s, v] for g, p, e, s, v in \
                zip(genes, probes, exp, sig, pval)]

    return Table(header=header, rows=rows)