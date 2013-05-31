import sys
sys.path.extend(['..'])

from cogent.util.unit_test import TestCase, main
import logging
from chippy.util.run_record import RunRecord, LOG_FN

class TestRunRecord(TestCase):
    def test_add_commands(self):
        """ test that RunRecord.addCommands correctly logs a max of three
            commands per line and doesn't fail with numbers that aren't
            multiples of three
        """
        logging.disable(logging.NOTSET)
        rr = RunRecord('test_add_commands')
        rr.addCommands([])
        rr.addCommands(['abc'])
        rr.addCommands(['abc', 'def'])
        rr.addCommands(['abc', 'def', 'ghi'])
        rr.addCommands(['abc', 'def', 'ghi', 'jkl'])
        rr.addCommands(['abc', 'def', 'ghi', 'jkl', 'mno'])

        recorded_lines = [
            'ChipPy.test_add_commands\tINFO\tcommand-line\tNo arguments given',
            'ChipPy.test_add_commands\tINFO\tabc command-line\tabc',
            'ChipPy.test_add_commands\tINFO\tabc command-line\tabc def',
            'ChipPy.test_add_commands\tINFO\tabc command-line\tabc def ghi',
            'ChipPy.test_add_commands\tINFO\tabc command-line\tabc def ghi',
            'ChipPy.test_add_commands\tINFO\tabc command-line\tjkl',
            'ChipPy.test_add_commands\tINFO\tabc command-line\tabc def ghi',
            'ChipPy.test_add_commands\tINFO\tabc command-line\tjkl mno'
        ]

        log_file = open(LOG_FN, 'r')
        for n, line in enumerate(log_file):
            line_parts = [lp.strip() for lp in line.split('\t')]
            #print repr(recorded_lines[n])
            #print repr('\t'.join(line_parts[1:]))
            assert '\t'.join(line_parts[1:]) == recorded_lines[n]

        logging.disable(logging.CRITICAL)

if __name__ == '__main__':
    main()