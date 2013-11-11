"""
    - probeset_label: label of the column containing probset id
    - ensembl_id_label: label of the column containing Ensembl Stable IDs
    - expression_label: label of the column containing absolute measure
        of expression
    - prob_label: label of the column containing raw probability of
        difference in gene expression between samples A/B
    - sig_label: label of the column classifying probabilities as
        significant after multiple test correction. 1 means up in A relative
        to B, -1 means down in A relative to B, 0 means no difference.
"""

__author__ = 'Cameron Jack'
__copyright__ = 'Copyright 2011-2013, Gavin Huttley, Cameron Jack, Anuj Pahwa'
__credits__ = ['Cameron Jack']
__license__ = 'GPL'
__maintainer__ = 'Cameron Jack'
__email__ = 'cameron.jack@anu.edu.au'
__status__ = 'pre-release'
__version__ = '0.1'

EXPR_HEADER = ['gene', 'probeset', 'exp']
DIFF_HEADER = EXPR_HEADER + ['sig', 'p_val']

