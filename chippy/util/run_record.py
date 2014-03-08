from cogent import LoadTable
__author__ = 'Cameron Jack, Gavin Huttley'
__copyright__ = 'Copyright 2011-2013, Gavin Huttley, Cameron Jack, Anuj Pahwa'
__credits__ = ['Gavin Huttley', 'Cameron Jack']
__license__ = 'GPL'
__maintainer__ = 'Cameron Jack'
__email__ = 'cameron.jack@anu.edu.au'
__status__ = 'Pre-release'
__version__ = '0.2'
import logging, sys, os, ConfigParser

LOG_FN = 'ChipPy.log'

class RunRecord(object):
    """ Python logging with extra, desired behaviours """
    def __init__(self, name=None):
        """
            Logger setup code from Python Logging Cookbook
            http://docs.python.org/2/howto/logging-cookbook.htm
        """
        super(RunRecord, self).__init__()

        base_name = 'ChipPy'
        if name is not None:
            self.logger = logging.getLogger(base_name+'.'+name)
            self.logger.propagate = False
        else:
            self.logger = logging.getLogger(base_name)

        # check for a log directory config file entry
        cfg = ConfigParser.ConfigParser()
        cfg.read('chippy.ini')
        try:
            home_dir = cfg.get('Log', 'Directory')
        except ConfigParser.NoSectionError, ConfigParser.NoOptionError:
            home_dir = None

        if home_dir is not None:
            self.log_path = os.path.join(home_dir, LOG_FN)
        else:
            self.log_path = LOG_FN

        if not len(self.logger.handlers):
            self.logger.setLevel(logging.DEBUG)
            # create file handler which logs even debug messages
            fh = logging.FileHandler(self.log_path)
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
            value_str = str(value)
            if len(value_str) == '':
                value_str = '<EMPTY>'
            self.logger.debug(message + '\t' + value_str)

    def addInfo(self, message='', value=''):
        """ adds a communication to the user about runtime """
        if type(value)==list or type(value)==tuple:
            for short_str in self._as_short_strings(value):
                self.logger.info(message + '\t' + short_str)
        else:
            value_str = str(value)
            if len(value_str) == '':
                value_str = '<EMPTY>'
            self.logger.info(message + '\t' + value_str)

    def addWarning(self, message='', value=''):
        """ adds a runtime warning """
        if type(value)==list or type(value)==tuple:
            for short_str in self._as_short_strings(value):
                self.logger.warning(message + '\t' + short_str)
        else:
            value_str = str(value)
            if len(value_str) == '':
                value_str = '<EMPTY>'
            self.logger.warning(message + '\t' + value_str)

    def addError(self, message='', value=''):
        """ report a non-critical runtime error """
        if type(value)==list or type(value)==tuple:
            for short_str in self._as_short_strings(value):
                self.logger.error(message + '\t' + short_str)
        else:
            value_str = str(value)
            if len(value_str) == '':
                value_str = '<EMPTY>'
            self.logger.error(message + '\t' + value_str)

    def addCritical(self, message='', value=''):
        """ report a catastrophic failure """
        if type(value)==list or type(value)==tuple:
            for short_str in self._as_short_strings(value):
                self.logger.critical(message + '\t' + short_str)
        else:
            value_str = str(value)
            if len(value_str) == '':
                value_str = '<EMPTY>'
            self.logger.critical(message + '\t' + value_str)

    def addCommands(self, command_args):
        """ groups together command-line args and adds to log """
        if len(command_args) == 0:
            self.addInfo('command-line', 'No arguments given')
        else:
            self.addInfo('command-line', command_args)

    def getMessageTable(self, last_n_lines=None, include_date=False):
        """
            Read the ChipPy.log file return as table, returning
            only the last n lines if passed an int.
        """

        log_file = open(self.log_path)
        records = []
        for line in log_file:
            line = line.strip()
            if len(line) > 0:
                if include_date:
                    records.append(line.split('\t')[0:])
                else:
                    records.append(line.split('\t')[1:]) # don't display date
        log_file.close()

        if records == []:
            return None

        if include_date:
            header = ['Date/time', 'code_block', 'level', 'message', 'value']
        else:
            header = ['code_block', 'level', 'message', 'value']

        if type(last_n_lines) is int: # return only last n lines of log file
            try:
                table = LoadTable(header=header, rows=records[-last_n_lines:], sep='\t')
            except IndexError:
                table = None
        else:
            try:
                table = LoadTable(header=header, rows=records, sep='\t')
            except IndexError:
                table = None
        return table

    def display(self):
        """ Displays the ChipPy.log file contents """
        table = self.getMessageTable(last_n_lines=30)
        if table is not None:
            print LOG_FN +' contains date-time stamps'
            print table

    def dieOnCritical(self, message='', value=''):
        """ log a CRITICAL failure, display all run records and die """
        self.addCritical(message=message, value=value)
        self.display()
        sys.exit(1)

    def empty_log(self):
        """ delete the contents of the log file """
        try:
            os.remove(self.log_path)
        except IOError:
            pass

        # recreate empty log file
        with open(self.log_path, 'w') as f:
            f.write('')
            f.close()

def remove_RR_log():
    """ delete the current ChipPy.log file, required when running tests """
    os.remove(LOG_FN)