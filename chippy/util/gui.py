import sys
from PyQt4 import QtGui

from argobs import ArgOb, DirPath, OpenFilePath, SaveFilePath,\
        ImportantStr, ImportantChoice, ImportantInt, ImportantFloat

__author__ = 'Cameron Jack'
__copyright__ = 'Copyright 2014, Cameron Jack'
__credits__ = ['Cameron Jack']
__license__ = 'GPL'
__version__ = '0.9.dev'
__maintainer__ = 'Cameron Jack'
__email__ = 'cameron.jack@anu.edu.au'
__status__ = 'Pre-release'

# TODO: have onLoad() populate the options
# TODO: create the command string from the options
# TODO: create mouse-over help strings and attach to widgets
# TODO: make the options window scrollable
# TODO: add in dividing lines for required and non-required options
# TODO: make sure more options support nargs!
# TODO: needs to validate user editable fields
# Note: There is no support for nested parsers or mutually exclusive options

class HelpText(object):
    """ Construction of argument information to text help for display """
    def __init__(self, long_form, num_args, short_form = '',
            type_or_action='', is_required='Not required'):
        self.num_args = num_args
        self.type_or_action = type_or_action
        self.is_required = is_required
        self.short_form = short_form
        self.long_form = long_form

    def constructHelpText(self):
        """ based on available info, join fields to form sentences """
        parts = [
            self.short_form,
            self.long_form,
            self.num_args,
            self.type_or_action,
            self.is_required
        ]
        return ' '.join(parts)

class AutoGUI(QtGui.QDialog):
    """
        Interprets ArgObs CLI argument objects and builds a GUI for the user
        and then returns a command-line for parsing.

        The basis of AutoGUI is ArgparseUi:
        https://pypi.python.org/pypi/argparseui
        but there is no requirement to use any particular parser and
        information for layouts and arguments comes from a purpose-built
        intermediate object layer named ArgObs.

        AutoGUI uses GridLayout to allow more flexibility with
        widget types and interface design - in particular the use of
        buttons to launch file and directory browser dialog boxes.
    """

    def __init__(self, argobs, window_title='Make your choice',
            header_text=None, footer_text=None, parent=None):
        """
            AutoGUI would generally be called from wherever CLI arguments
            are defined as ArgObs.
        """
        super(AutoGUI, self).__init__(parent)
        self.setWindowTitle(window_title)
        self.argobs = argobs
        self.header_text = header_text
        self.footer_text = footer_text

        self.config_fn = None

        self.mainLayout = QtGui.QVBoxLayout(self)
        self.setLayout(self.mainLayout)
        self.createContents()

    def createContents(self):
        """ Fill the main layout with joy and stuff """
        self.addHeaderToGUI()
        self.createOptions()
        self.addFooterToGUI()
        self.addButtonsToGUI()

    def addHeaderToGUI(self):
        """
            Set descriptive text above the options.
            Equivalent to argparse.description
        """
        self.header = QtGui.QWidget(self)
        self.headerLayout = QtGui.QVBoxLayout(self.header)
        self.header.setLayout(self.headerLayout)

        if self.header_text is None:
            self.header_text = '<b>Choose which options to include, '+\
                    'and define their value below</b>'
        label = QtGui.QLabel(self.header_text)
        self.headerLayout.addWidget(label)
        # Now add to mainLayout
        self.mainLayout.addWidget(self.header)

    def makeSelectDirEvent(self, line_edit):
        def selectDirEvent():
            line_edit.setText(QtGui.QFileDialog.getExistingDirectory())
        return selectDirEvent

    def makeSelectOpenFileEvent(self, line_edit, nargs=None, filter_=''):
        def selectOpenFileEvent():
            if nargs is None:
                line_edit.setText(QtGui.QFileDialog.\
                        getOpenFileName(filter=filter_))
            else:
                line_edit.addItems(QtGui.QFileDialog.\
                        getOpenFileNames(filter=filter_))
        return selectOpenFileEvent

    def makeSelectSaveFileEvent(self, line_edit, filter_=''):
        def selectSaveFileEvent():
            line_edit.setText(QtGui.QFileDialog.\
                    getSaveFileName(filter=filter_))
        return selectSaveFileEvent

    def _buildFileButtonLineEdit(self, arg):
        """
            Add a button widget (which invokes a file browser dialog)
            and a lineEdit widget to gridlayout.
            Used for all derivatives of ArgOb.FilePath
        """
        button = QtGui.QPushButton(arg.long_form.lstrip('-'))
        if arg.nargs is None:
            line_edit = QtGui.QLineEdit()
        else:
            line_edit = QtGui.QListWidget()
        if arg.arg_type is OpenFilePath:
            button.clicked.connect(self.makeSelectOpenFileEvent(line_edit,
                    nargs=arg.nargs, filter_=arg.filter_))
        elif arg.arg_type is SaveFilePath:
            button.clicked.connect(self.makeSelectSaveFileEvent(line_edit))
        elif arg.arg_type is DirPath:
            button.clicked.connect(self.makeSelectDirEvent(line_edit))
        else:
            QtGui.QMessageBox.critical(self, 'Critical',
                'Invalid arg_type {0}'.format(str(arg.arg_type)))

        self.optionsLayout.addWidget(button, self.current_row, 0)
        self.optionsLayout.addWidget(line_edit, self.current_row, 1)
        self.current_row += 1

    def _buildLabelCheckbox(self, arg):
        """
            Add a label widget and a checkbox widget to gridlayout.
            Used for action.
        """
        checkbox = QtGui.QCheckBox()
        if arg.default or arg.action == 'store_false':
            checkbox.setChecked(True)
        label = QtGui.QLabel(arg.long_form.lstrip('-'))
        #label.setBuddy(checkbox)

        self.optionsLayout.addWidget(label, self.current_row, 0)
        self.optionsLayout.addWidget(checkbox, self.current_row, 1)
        self.current_row += 1

    def _buildLabelCombobox(self, arg):
        """
            Add a label widget and a combobox widget to gridlayout.
            Used for choice, ImportantChoice
        """
        combobox = QtGui.QComboBox()
        label = QtGui.QLabel(arg.long_form.lstrip('-'))
        #label.setBuddy(combobox)
        for c in arg.choices:
            combobox.addItem('{0}'.format(c))
        if arg.default is not None:
            combobox.setCurrentIndex(combobox.findText('{0}'.\
                    format(arg.default)))

        self.optionsLayout.addWidget(label, self.current_row, 0)
        self.optionsLayout.addWidget(combobox, self.current_row, 1)
        self.current_row += 1

    def _buildLabelLineEdit(self, arg):
        """
            Add a label widget and a lineEdit widget to gridlayout.
            Used for str, numbers, ImportantStr, ImportantNumbers
        """
        line_edit = QtGui.QLineEdit()
        label = QtGui.QLabel(arg.long_form.lstrip('-'))
        #label.setBuddy(line_edit)

        self.optionsLayout.addWidget(label, self.current_row, 0)
        self.optionsLayout.addWidget(line_edit, self.current_row, 1)
        self.current_row += 1

    def _createOrderedWidgets(self, argobs):
        """
            Build widgets appropriate to the arg type, and then add them to
            the interface. Widgets should be displayed in order by the
            arg_type: DirPath should come first, then FilePath,
            ImportantStr, actions, choices, numbers and finally strings.
        """
        row = 0 # gridlayout row
        # group argobs by type so they are placed together
        for arg in [a for a in argobs if a.arg_type is DirPath]:
            self._buildFileButtonLineEdit(arg)
        for arg in [a for a in argobs if a.arg_type is OpenFilePath]:
            self._buildFileButtonLineEdit(arg)
        for arg in [a for a in argobs if a.arg_type is SaveFilePath]:
            self._buildFileButtonLineEdit(arg)
        for arg in [a for a in argobs if a.arg_type is ImportantChoice]:
            self._buildLabelCombobox(arg)
        for arg in [a for a in argobs if a.arg_type is ImportantInt]:
            self._buildLabelLineEdit(arg)
        for arg in [a for a in argobs if a.arg_type is ImportantFloat]:
            self._buildLabelLineEdit(arg)
        for arg in [a for a in argobs if a.arg_type is ImportantStr]:
            self._buildLabelLineEdit(arg)
        for arg in [a for a in argobs if a.action is not None]:
            self._buildLabelCheckbox(arg)
        for arg in [a for a in argobs if a.choices is not None]:
            self._buildLabelCombobox(arg)
        for arg in [a for a in argobs if a.arg_type is int]:
            self._buildLabelLineEdit(arg)
        for arg in [a for a in argobs if a.arg_type is float]:
            self._buildLabelLineEdit(arg)
        for arg in [a for a in argobs if a.arg_type is str]:
            self._buildLabelLineEdit(arg)

    def createOptions(self):
        """
            Create option names (long form only) and their input field
            as appropriate for their type and action. Files should open
            a browser dialog so that the user can easily select the
            locations. Full argument information should be available as
            mouse-over help.
        """
        self.options = QtGui.QWidget(self)
        self.optionsLayout = QtGui.QGridLayout(self.options)
        self.options.setLayout(self.optionsLayout)

        # Split args into required and non-required options
        required_args = []
        nonreq_args = []
        for arg in self.argobs:
            if arg.display:
                if arg.required:
                    required_args.append(arg)
                else:
                    nonreq_args.append(arg)
        self.current_row = 0
        # need to add a "Required inputs:" line
        self._createOrderedWidgets(required_args)
        # need to add a "Optional inputs:" line
        self._createOrderedWidgets(nonreq_args)

        self.mainLayout.addWidget(self.options)

    def addButton(self, label, do_func):
        """ creates an individual button, that does func """
        button = QtGui.QPushButton(label, self.buttons)
        button.clicked.connect(do_func)
        self.buttonsLayout.addWidget(button)
        self.buttonsLayout.addSpacerItem(QtGui.QSpacerItem(\
                20, 1, QtGui.QSizePolicy.Expanding))
        return button

    def addButtonsToGUI(self):
        """
            Create Load, Ok, Save, Save As, Reset and Cancel option buttons
        """
        self.buttons = QtGui.QWidget(self)
        self.buttonsLayout = QtGui.QHBoxLayout(self.buttons)
        self.buttons.setLayout(self.buttonsLayout)

        self.ok_button = self.addButton('Ok', self.onOk)
        self.load_button = self.addButton('Load', self.onLoad)
        self.save_button = self.addButton('Save', self.onSave)
        self.save_as_button = self.addButton('Save as', self.onSaveAs)
        self.reset_button = self.addButton('Reset', self.onReset)
        self.cancel_button = self.addButton('Cancel', self.onCancel)

        self.mainLayout.addWidget(self.buttons)

    def addFooterToGUI(self):
        """
            Set any footer text below the options and above the buttons.
            Equivalent to argparse.epilog
        """
        self.footer = QtGui.QWidget(self)
        self.footerLayout = QtGui.QVBoxLayout(self.footer)
        self.footer.setLayout(self.footerLayout)

        if self.footer_text is None:
            self.footer_txt = 'This options dialog is auto-generated by autoGUI'
        label = QtGui.QLabel(self.footer_txt)
        self.footerLayout.addWidget(label)
        # Now add to mainLayout
        self.mainLayout.addWidget(self.footer)

    def createCommandLine(self):
        """ glue together the arg.long_names and the user inputs """
        pass

    # Button press functions

    def onOk(self):
        """ create command-line from form and return it """
        self.accept()

    def onLoad(self):
        """ Create a dialog to open a saved option configuration """
        pass

    def onSave(self):
        """ Save current option configuration to file """
        if not self.config_fn:
            self.onSaveAs()
        else:
            try:
                with open(self.config_fn, 'w') as f:
                    c = self.makeCommandLine()
                    for line in c:
                        f.write(line+'\n')
            except IOError:
                import os.path
                QtGui.QMessageBox.critical(self, 'Critical',
                        "Couldn't write to file {0}".format\
                        (os.path.abspath(self.filename)))

    def onSaveAs(self):
        """ Create a dialog to save option configuration """
        config_fn = QtGui.QFileDialog.getSaveFileName()
        if config_fn:
            print config_fn
            self.config_fn = config_fn
            self.onSave()

    def onReset(self):
        """ Clear all current option fields """
        response = QtGui.QMessageBox.question(self, 'Alert',
                'Are you sure you want to clear all inputs?',
                QtGui.QMessageBox.Yes | QtGui.QMessageBox.No,
                QtGui.QMessageBox.No )
        if response == QtGui.QMessageBox.Yes:
            if self.mainLayout is not None:
                while self.mainLayout.count():
                    item = self.mainLayout.takeAt(0)
                    widget = item.widget()
                    if widget is not None:
                        widget.deleteLater()
                    else:
                        self.clearLayout(item.layout())
            self.createContents()

    def onCancel(self):
        """ Close the application and return NoneType """
        self.reject()

def main():
    app = QtGui.QApplication(sys.argv)
    args = [
        ArgOb('--num_genes', type=int,
                help='Number of ranked genes to use in study'),
        ArgOb('--group_location', default='all',
                choices=['all', 'top', 'middle', 'bottom'],
                help='The representative group in a study to form a plot line'),
        ArgOb('--ranks', action='store_true',
                help='Use rank-based data values instead of counts or '+\
                'expression'),
        ArgOb('--sample_extremes', type=float,
                default=0.0, help='Proportion of least and most absolute '+\
                'expressed genes to treat separately. Set to 0.0 to disable.'),
        ArgOb('--name', help='some name', required=True),
        ArgOb('--in_file', type=OpenFilePath, nargs='+', help='path to infile'),
        ArgOb('--some_dir', type=DirPath, help='path to directory'),
        ArgOb('--save_file', type=SaveFilePath, help='path save file')
    ]

    a = AutoGUI(args)
    a.show()
    app.exec_()
    print ("Ok" if a.result() == 1 else "Cancel")
    if a.result() == 1: # Ok pressed
        return a.parse_args()
    else:
        print 'Execution cancelled.'
        sys.exit(0)

if __name__ == '__main__':
    main()