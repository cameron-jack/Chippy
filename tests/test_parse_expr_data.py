import sys
sys.path.extend(['..'])

from cogent import LoadTable
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files

from chippy.parse.expr_data import _check_expr_headers, _check_diff_headers,\
        _validate_probes_scores, _remove_multimapped_probesets

_sample_dump = LoadTable(header=['ENSEMBL', 'probeset', 'exp'],
        rows=[['id1',"0|1|2","13.6|13.4|13.6"],
        ['id2',"3|1","9.9|13.6"], # this gene should be lost when filtered
        ['id3',"4|5","12.7|13.4"],
        ['id4',"6","13.4"],
        ['id5',"7|8|3","6.0|6.0|4.5"],
        ['id6',"9|10|11|12","5.4|6.8|6.6|6.2"],
        ['id8',"13","12.7"],
        ['id9',"14","12.7"],
        ['id10',"15","12.7"]])

class ExcludingProbesets(TestCase):
    """test that excluding probesets works correctly"""

    def test_check_expr_headers(self):
        """ check that headers are identified corrects, as is the presence/
            absence of a probeset column label. Make sure it fails if the
            columns are incorrectly ordered or labelled.
        """
        header_row = ['ENSEMBL', 'probeset', 'exp']
        gene_col, probe_col, exp_col, probes_present = _check_expr_headers(
                header_row, stable_id_label='ENSEMBL',
                probeset_label='probeset', exp_label='exp')
        self.assertTrue(probes_present)
        self.assertEqual(gene_col, 0)
        self.assertEqual(probe_col, 1)
        self.assertEqual(exp_col, 2)

        header_row = ['ENSEMBL', 'exp']
        gene_col, probe_col, exp_col, probes_present = _check_expr_headers(
                header_row, stable_id_label='ENSEMBL',
                probeset_label='probeset', exp_label='exp')
        self.assertFalse(probes_present)
        self.assertEqual(gene_col, 0)
        self.assertEqual(exp_col, 1)

        # will still work but think probeset is missing
        header_row = ['ENSEMBL', 'exp', 'probeset']
        gene_col, probe_col, exp_col, probes_present = _check_expr_headers(
                header_row, stable_id_label='ENSEMBL',
                probeset_label='probeset', exp_label='exp')
        self.assertTrue(probes_present)
        self.assertEqual(gene_col, 0)
        self.assertEqual(probe_col, 2)
        self.assertEqual(exp_col, 1)

        # will fail because gene label is wrong
        header_row = ['gene', 'probeset', 'exp']
        self.assertRaises(SystemExit, _check_expr_headers, header_row,
                stable_id_label='ENSEMBL', probeset_label='probeset',
                exp_label='exp')

        # will not because probeset label is wrong
        header_row = ['ENSEMBL', 'probes', 'exp']
        gene_col, probe_col, exp_col, probes_present = _check_expr_headers(
                header_row, stable_id_label='ENSEMBL',
                probeset_label='probeset', exp_label='exp')
        self.assertFalse(probes_present)
        self.assertEqual(gene_col, 0)
        self.assertEqual(exp_col, 2)

        # will fail because exp label is wrong
        header_row = ['ENSEMBL', 'probeset', 'scores']
        self.assertRaises(SystemExit, _check_expr_headers, header_row,
                stable_id_label='ENSEMBL', probeset_label='probeset',
                exp_label='exp')

    def test_check_diff_headers(self):
        """ check the extra column identifiers in a diff file, with/without
            the presence of probesets. True/False are for the 'probes_present'
            arg.
        """
        header_row = ['gene', 'probeset', 'exp', 'sig', 'rawp']
        sig_col, pval_col = _check_diff_headers(header_row, sig_label='sig',
                pval_label='rawp')
        self.assertEqual(sig_col, 3)
        self.assertEqual(pval_col, 4)

        header_row = ['gene', 'exp', 'sig', 'rawp']
        sig_col, pval_col = _check_diff_headers(header_row, sig_label='sig',
                pval_label='rawp')
        self.assertEqual(sig_col, 2)
        self.assertEqual(pval_col, 3)

        header_row = ['gene', 'exp', 'sig', 'rawp']
        self.assertRaises(SystemExit, _check_diff_headers, header_row,
                sig_label='signif', pval_label='rawp')

        header_row = ['gene', 'probeset', 'exp', 'sig', 'rawp']
        self.assertRaises(SystemExit, _check_diff_headers, header_row,
                sig_label='signif', pval_label='p_vals')

    def test_validate_probes_scores(self):
        """ make sure that probes/scores are matched and if not then the
            whole gene entry is nuked.
        """

        genes = ['id1', 'id2', 'id3']
        probes = [['0','1','2'], ['3','1'], ['4','5']]
        scores = [[13.6, 13.4,13.6], [9.9, 13.6], [12.7, 13.4]]
        v_genes, v_probes, v_scores = _validate_probes_scores(genes,
                probes, scores)
        self.assertEqual(genes, v_genes)
        self.assertEqual(probes, v_probes)
        self.assertEqual(scores, v_scores)

        genes = ['id1', 'id2', 'id3']
        probes = [['0','1','2'], ['3'], ['4','5']]
        scores = [[13.6, 13.4,13.6], [9.9, 13.6], [12.7]]
        v_genes, v_probes, v_scores = _validate_probes_scores(genes,
            probes, scores)

        self.assertEqual(v_genes, [genes[0]])
        self.assertEqual(v_probes, [probes[0]])
        self.assertEqual(v_scores, [scores[0]])

    def test_remove_multimapped_probesets(self):
        """
            If the same probeset shows up in more than 1 gene then the
            probeset and its matching score need to be removed from every
            gene that shares it.
        """
        genes = ['id1', 'id2', 'id3', 'id4', 'id5']
        probes = [['0','1','2'], ['3','1'], ['4','5'], ['1','3'], ['0']]
        scores = [[13.6, 13.4,13.6], [9.9, 13.6], [12.7, 13.4], [10.3, 11.7],
                [12.1]]

        v_genes, v_probes, v_scores = _remove_multimapped_probesets(genes,
                probes, scores)

        expected_genes = ['id1', 'id3']
        expected_probes = [['2'], ['4','5']]
        expected_scores = [[13.6], [12.7, 13.4]]

        self.assertEqual(v_genes, expected_genes)
        self.assertEqual(v_probes, expected_probes)
        self.assertEqual(v_scores, expected_scores)

if __name__ == '__main__':
    main()

