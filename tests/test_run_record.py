import sys
sys.path.extend(['..'])

from cogent.util.unit_test import TestCase, main
import logging
from chippy.util.run_record import RunRecord, LOG_FN

class TestRunRecord(TestCase):
    def test_add_commands(self):
        """ test that RunRecord.addCommands correctly logs long lines of text
        """
        logging.disable(logging.NOTSET)
        rr = RunRecord('test_add_commands')
        rr.addCommands([])
        cmd_line = 'This is a list of command arguments that probably '+\
                   'do not exist in the real world'
        cmds = cmd_line.split(' ')
        rr.addCommands(cmds)

        recorded_lines = [
            'ChipPy.test_add_commands\tINFO\tcommand-line\tNo arguments given',
            'ChipPy.test_add_commands\tINFO\tcommand-line\tThis is a list of command arguments',
            'ChipPy.test_add_commands\tINFO\tcommand-line\tthat probably do not exist in the real',
            'ChipPy.test_add_commands\tINFO\tcommand-line\tworld'
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