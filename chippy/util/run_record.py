from cogent import LoadTable
__author__ = "Gavin Huttley, Cameron Jack"
__copyright__ = "Copyright 2012, Gavin Huttley, Cameron Jack, Anuj Pahwa"
__credits__ = ["Gavin Huttley", "Cameron Jack"]
__license__ = "GPL"
__maintainer__ = "Cameron Jack"
__email__ = "cameron.jack@anu.edu.au"
__status__ = "Pre-release"
__version__ = '0.2'
import logging, sys, os
from math import ceil

LOG_FN = 'ChipPy.log'

class RunRecord(object):
    """ Python logging with legacy behaviours """
    def __init__(self, name=None):
        super(RunRecord, self).__init__()

        # Logger setup code from Python Logging Cookbook
        # http://docs.python.org/2/howto/logging-cookbook.html
        base_name = 'ChipPy'
        if name is not None:
            self.logger = logging.getLogger(base_name+'.'+name)
        else:
            self.logger = logging.getLogger(base_name)
        if not len(self.logger.handlers):
            self.logger.setLevel(logging.DEBUG)
            # create file handler which logs even debug messages
            fh = logging.FileHandler(LOG_FN)
            fh.setLevel(logging.DEBUG)
            # create console handler with a higher log level
            ch = logging.StreamHandler(sys.stderr)
            ch.setLevel(logging.ERROR)
            # create formatter and add it to the handlers
            formatter = logging.Formatter('%(asctime)s \t %(name)s \t %(levelname)s \t %(message)s')
            fh.setFormatter(formatter)
            formatter = logging.Formatter('%(name)s \t %(levelname)s \t %(message)s')
            ch.setFormatter(formatter)
            # add the handlers to the logger
            self.logger.addHandler(fh)
            self.logger.addHandler(ch)

    def _as_short_strings(self, value_parts):
        """ splits a multi-component value (e.g. list) into
            parts of no more than max_chars
        """
        MAX_CHARS = 35
        print_parts = ''
        for part in value_parts:
            if len(print_parts) < MAX_CHARS:
                print_parts = ' '.join([print_parts, str(part)])
            else:
                yield print_parts
                print_parts = str(part)
        yield print_parts

    def addDebug(self, message='', value=''):
        """ adds a debugging message about an execution """
        if type(value)==list or type(value)==tuple:
            for short_str in self._as_short_strings(value):
                self.logger.debug(message + '\t' + short_str)
        else:
            self.logger.debug(message + '\t' + str(value))

    def addInfo(self, message='', value=''):
        """ adds a communication to the user about runtime """
        if type(value)==list or type(value)==tuple:
            for short_str in self._as_short_strings(value):
                self.logger.info(message + '\t' + short_str)
        else:
            self.logger.info(message + '\t' + str(value))

    def addWarning(self, message='', value=''):
        """ adds a runtime warning """
        if type(value)==list or type(value)==tuple:
            for short_str in self._as_short_strings(value):
                self.logger.warning(message + '\t' + short_str)
        else:
            self.logger.warning(message + '\t' + str(value))

    def addError(self, message='', value=''):
        """ report a non-critical runtime error """
        if type(value)==list or type(value)==tuple:
            for short_str in self._as_short_strings(value):
                self.logger.error(message + '\t' + short_str)
        else:
            self.logger.error(message + '\t' + str(value))

    def addCritical(self, message='', value=''):
        """ report a catastrophic failure """
        if type(value)==list or type(value)==tuple:
            for short_str in self._as_short_strings(value):
                self.logger.critical(message + '\t' + short_str)
        else:
            self.logger.critical(message + '\t' + str(value))

    def addCommands(self, command_args):
        """ groups together command-line args and adds to log """
        if len(command_args) == 0:
            self.addInfo('command-line', 'No arguments given')
        else:
            self.addInfo('command-line', command_args)

    def getMessageTable(self):
        """docstring for display"""

        log_file = open(LOG_FN)
        records = []
        for line in log_file:
            line = line.strip()
            if len(line) > 0:
                records.append(line.split('\t')[1:]) # don't display date
        log_file.close()

        #header = ['Date/time', 'code_block', 'level', 'message', 'value']
        header = ['code_block', 'level', 'message', 'value']
        table = LoadTable(header=header, rows=records, sep='\t')
        return table

    def display(self):
        table = self.getMessageTable()
        if table is None:
            print 'No log messages'
        else:
            print LOG_FN +' contains date-time stamps'
            print table

    def dieOnCritical(self, message='', value=''):
        """ log a CRITICAL failure, display all run records and die """
        self.addCritical(message=message, value=value)
        self.display()
        sys.exit(1)

def remove_RR_log():
    """ delete the current ChipPy.log file, such as when running tests """
    os.remove(LOG_FN)