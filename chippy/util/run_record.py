from cogent import LoadTable
__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"
__version__ = '0.1'

class RunRecord(object):
    """object for recording program messages"""
    def __init__(self):
        super(RunRecord, self).__init__()
        self.records = []
        
    def addMessage(self, program_name, error_type, message, value):
        """add a message about an execution"""
        self.records.append([program_name, error_type, message, value])
    
    def addInfo(self, program_name, message, value):
        self.addMessage(program_name, 'info', message, value)
    
    def addWarning(self, program_name, message, value):
        self.addMessage(program_name, 'warning', message, value)
    
    def addError(self, program_name, message, value):
        self.addMessage(program_name, 'error', message, value)
    
    def getMessageTable(self):
        """docstring for display"""
        if not self.records:
            return None
        
        header = ['program_name', 'message type', 'message', 'value']
        table = LoadTable(header=header, rows=self.records)
        return table
    
    def display(self):
        table = self.getMessageTable()
        if table is None:
            print 'No log messages'
        else:
            print table
        
    

