import sys
from PyQt4 import QtGui

from argobs import ArgOb, DirPath, OpenFilePath, SaveFilePath

__author__ = 'Cameron Jack'
__copyright__ = 'Copyright 2014, Cameron Jack'
__credits__ = ['Cameron Jack']
__license__ = 'GPL'
__version__ = '0.9.dev'
__maintainer__ = 'Cameron Jack'
__email__ = 'cameron.jack@anu.edu.au'
__status__ = 'Pre-release'

# Note: There is no support for nested parsers or mutually exclusive options
# Not all types support nargs > 1
# Configuration load/save is not properly implemented and available

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

        # Bring this GUI to the front
        self.show()
        self.raise_()
        self.activateWindow()

    def createContents(self):
        """ Fill the main layout with a scrollable widget """

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

    def _buildLabelLineEditFileButton(self, arg, layout):
        """
            Add a button widget (which invokes a file browser dialog)
            and a lineEdit widget to gridlayout.
            Used for all derivatives of ArgOb.FilePath
        """
        label = QtGui.QLabel(arg.long_form.lstrip('-').replace('_', ' '))
        button = QtGui.QToolButton()
        if sys.platform.startswith('darwin'):
            button.setText('Browse')
        else:
            button.setIcon(QtGui.QIcon().fromTheme('folder-new'))

        #button.setAutoDefault(False)
        if arg.nargs is None:
            line_edit = QtGui.QLineEdit()
        else:
            line_edit = QtGui.QListWidget()
        if arg.arg_type is OpenFilePath:
            button.clicked.connect(self.makeSelectOpenFileEvent(line_edit,
                    nargs=arg.nargs, filter_=arg.filter_))
            button.setToolTip('Browse for file(s) to open')
        elif arg.arg_type is SaveFilePath:
            button.clicked.connect(self.makeSelectSaveFileEvent(line_edit))
            button.setToolTip('Browse for save file path')
        elif arg.arg_type is DirPath:
            button.clicked.connect(self.makeSelectDirEvent(line_edit))
            button.setToolTip('Browse for directory to select')
        else:
            QtGui.QMessageBox.critical(self, 'Critical',
                    'Invalid arg_type {0}'.format(str(arg.arg_type)))
        self._addToLayout(arg, layout, [label, line_edit, button])

    def _buildLabelCheckbox(self, arg, layout):
        """
            Add a label widget and a checkbox widget to gridlayout.
            Used for action.
        """
        #checkbox = QtGui.QCheckBox()
        #if arg.default or arg.action == 'store_false':
        #    checkbox.setChecked(True)
        label = QtGui.QLabel(arg.long_form.lstrip('-').replace('_', ' '))
        self._addToLayout(arg, layout, [label]) #, checkbox])

    def _buildLabelCombobox(self, arg, layout):
        """
            Add a label widget and a combobox widget to gridlayout.
            Used for choice, ImportantChoice
        """
        combobox = QtGui.QComboBox()
        label = QtGui.QLabel(arg.long_form.lstrip('-').replace('_', ' '))
        for c in arg.choices:
            combobox.addItem('{0}'.format(c))
        if arg.default is not None:
            combobox.setCurrentIndex(combobox.findText('{0}'.\
                    format(arg.default)))
        self._addToLayout(arg, layout, [label, combobox])

    def _buildLabelDoubleSpinBox(self, arg, layout):
        """
            Add a label widget and a DoubleSpinBox widget to gridlayout.
            Used for float, ImportantFloat.
        """
        double_spinbox = QtGui.QDoubleSpinBox()
        double_spinbox.setDecimals(10)
        double_spinbox.setRange(-999999999.0,999999999.0)
        if arg.default is not None:
            double_spinbox.setValue(arg.default)
        label = QtGui.QLabel(arg.long_form.lstrip('-').replace('_', ' '))
        self._addToLayout(arg, layout, [label, double_spinbox])

    def _buildLabelSpinBox(self, arg, layout):
        """
            Add a label widget and a SpinBox widget to gridlayout.
            Used for str, numbers, ImportantStr, ImportantNumbers
        """
        spinbox = QtGui.QSpinBox()
        spinbox.setRange(-999999999,999999999)
        if arg.default is not None:
            spinbox.setValue(arg.default)
        label = QtGui.QLabel(arg.long_form.lstrip('-').replace('_', ' '))
        self._addToLayout(arg, layout, [label, spinbox])

    def _buildLabelLineEdit(self, arg, layout):
        """
            Add a label widget and a lineEdit widget to gridlayout.
            Used for str, numbers, ImportantStr, ImportantNumbers
        """
        label = QtGui.QLabel(arg.long_form.lstrip('-').replace('_', ' '))
        line_edit = QtGui.QLineEdit()
        if arg.default is not None:
            line_edit.setText(arg.default)
        self._addToLayout(arg, layout, [label, line_edit])

    def _addToLayout(self, arg, layout, widgets):
        """ set up all widgets and add them to the layout """
        col = 0
        if not arg.required:
            include = QtGui.QCheckBox()
            layout.addWidget(include, self.current_row, col)
            col += 1
            for widget in widgets:
                if len(widget.toolTip()) == 0:
                    widget.setToolTip(arg.gui_help)
                self.disableOnClick(widget)(False)
                include.clicked.connect(self.disableOnClick(widget))
                layout.addWidget(widget, self.current_row, col)
                col += 1
        else:
            col = 1
            for widget in widgets:
                widget.setToolTip(arg.gui_help)
                layout.addWidget(widget, self.current_row, col)
                col += 1
        self.current_row += 1

    def disableOnClick(self, widget):
        """ Factory to disable/enable widgets """
        def disable(state):
            widget.setEnabled(state)
        return disable

    def _dispatchWidgetBuilding(self, argobs, layout):
        """ group argobs by type so they are placed together """
        for arg in [a for a in argobs if a.arg_type is DirPath]:
            self._buildLabelLineEditFileButton(arg, layout)
        for arg in [a for a in argobs if a.arg_type is OpenFilePath]:
            self._buildLabelLineEditFileButton(arg, layout)
        for arg in [a for a in argobs if a.arg_type is SaveFilePath]:
            self._buildLabelLineEditFileButton(arg, layout)
        for arg in [a for a in argobs if a.action is not None]:
            self._buildLabelCheckbox(arg, layout)
        for arg in [a for a in argobs if a.choices is not None]:
            self._buildLabelCombobox(arg, layout)
        for arg in [a for a in argobs if a.arg_type is int]:
            self._buildLabelSpinBox(arg, layout)
        for arg in [a for a in argobs if a.arg_type is float]:
            self._buildLabelDoubleSpinBox(arg, layout)
        for arg in [a for a in argobs if a.arg_type is str]:
            self._buildLabelLineEdit(arg, layout)

    def _createOrderedWidgets(self, argobs, layout):
        """
            Build widgets appropriate to the arg type, and then add them to
            the interface. Widgets should be displayed in order by the
            arg_type: DirPath should come first, then FilePath,
            ImportantStr, actions, choices, numbers and finally strings.
        """
        for arg in argobs:
            arg.gui_help = self.createMouseOverHelp(arg)

        self._dispatchWidgetBuilding([a for a in argobs if a.important],
                layout)
        self._dispatchWidgetBuilding([a for a in argobs if not a.important],
                layout)

    def createOptions(self):
        """ Options are split between required and optional """
        # master, scrollable panel
        self.options = QtGui.QWidget(self)
        self.optionsLayout = QtGui.QVBoxLayout(self.options)
        self.options.setLayout(self.optionsLayout)

        self.required = QtGui.QGroupBox(title='Required', flat=True)
        self.required.setToolTip('Required fields')
        self.requiredLayout = QtGui.QGridLayout(self.required)
        self.required.setLayout(self.requiredLayout)

        self.optional = QtGui.QGroupBox(title='Optional', flat=True)
        self.optional.setToolTip('Optional fields')
        self.optionalLayout = QtGui.QGridLayout(self.optional)
        self.optional.setLayout(self.optionalLayout)

        self.optionsLayout.addWidget(self.required)
        self.optionsLayout.addWidget(self.optional)

        # Split args into required and non-required options
        required_args = []
        optional_args = []
        for arg in self.argobs:
            if arg.display:
                if arg.required:
                    required_args.append(arg)
                else:
                    optional_args.append(arg)
        self.current_row = 0
        self._createOrderedWidgets(required_args, self.requiredLayout)
        self.current_row = 0
        self._createOrderedWidgets(optional_args, self.optionalLayout)

        self.mainLayout.addWidget(self.options)

        self.scrollableArea = QtGui.QScrollArea(self)
        self.scrollableArea.setWidgetResizable(True)
        self.scrollableArea.setFrameShape(QtGui.QFrame.NoFrame)
        self.scrollableArea.setWidget(self.options)
        self.mainLayout.addWidget(self.scrollableArea)

    def addCoreButton(self, label, do_func):
        """ creates an individual button, that does func """
        button = QtGui.QPushButton(label, self.buttons)
        button.setAutoDefault(False)
        button.clicked.connect(do_func)
        self.buttonsLayout.addWidget(button)
        self.buttonsLayout.addSpacerItem(QtGui.QSpacerItem(\
                20, 1, QtGui.QSizePolicy.Expanding))
        return button

    def addButtonsToGUI(self):
        """
            Create Load, Ok, Save, Save As, Reset and Cancel option buttons
            Load, Save, Save As are probably unnecessary
        """
        self.buttons = QtGui.QWidget(self)
        self.buttonsLayout = QtGui.QHBoxLayout(self.buttons)
        self.buttons.setLayout(self.buttonsLayout)

        self.ok_button = self.addCoreButton('Ok', self.onOk)
        #self.load_button = self.addCoreButton('Load', self.onLoad)
        #self.save_button = self.addCoreButton('Save', self.onSave)
        #self.save_as_button = self.addCoreButton('Save as', self.onSaveAs)
        self.reset_button = self.addCoreButton('Reset', self.onReset)
        self.cancel_button = self.addCoreButton('Cancel', self.onCancel)
        self.cancel_button.setFocus()

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
        self.mainLayout.addWidget(self.footer)

    def _checkNameValue(self, label, value_widget):
        """
            Check label against arg.long_name and type check value_widget
        """
        label = label.widget()
        label = str(label.text()).replace(' ', '_')
        # Handle cases with no checkbox
        if value_widget is None:
            for arg in self.argobs:
                if label in arg.long_form:
                    return arg.long_form, None
        else:
            value_widget = value_widget.widget()

        if 'QCheckBox' in str(type(value_widget)):
            value = value_widget.checkState()
        elif 'QLineEdit' in str(type(value_widget)):
            value = value_widget.text()
        elif 'QComboBox' in str(type(value_widget)):
            value = value_widget.currentText()
        elif 'QSpinBox' in str(type(value_widget)):
            value = value_widget.value()
        elif 'QDoubleSpinBox' in str(type(value_widget)):
            value = value_widget.value()
        elif 'QListWidget' in str(type(value_widget)):
            items = []
            for i in xrange(value_widget.count()):
                items.append(value_widget.item(i))
            value = [str(i.text()) for i in items]
        else:
            print 'uh oh! New widget type needed!'
            print value_widget, type(value_widget)

        # now check against original type to ensure correct formatting
        name = None
        input = None
        for arg in self.argobs:
            if label in arg.long_form:
                name = arg.long_form
                # no need to add anything for Action
                if type(value) is list:
                    input = ', '.join(value)
                else:
                    if arg.arg_type is not None:
                        input = str(value)
                    elif arg.choices is not None:
                        input = str(value)
        return name, input

    def makeCommandLine(self):
        """ glue together the arg.long_names and the user inputs """
        optional_parts = []
        positional_parts = []

        if self.requiredLayout.rowCount() > 0:
            for row in range(self.requiredLayout.rowCount()):
                if self.requiredLayout.itemAtPosition(row,2) is not None:
                    name, input = self._checkNameValue(\
                            self.requiredLayout.itemAtPosition(row, 1),
                            self.requiredLayout.itemAtPosition(row, 2))
                    if name.startswith('-'):
                        optional_parts.append(name)
                        if input is not None:
                            optional_parts.append(input)
                    else: # positional args have no name
                        if input is not None:
                            positional_parts.append(input)

        if self.optionalLayout.rowCount() > 0:
            for row in range(self.optionalLayout.rowCount()):
                include = self.optionalLayout.itemAtPosition(row, 0).widget()
                if include.checkState():
                    name, input = self._checkNameValue(\
                        self.optionalLayout.itemAtPosition(row, 1),
                        self.optionalLayout.itemAtPosition(row, 2))
                    if name.startswith('-'):
                        optional_parts.append(name)
                        if input is not None:
                            optional_parts.append(input)
                    else: # positional args have no name
                        if input is not None:
                            positional_parts.append(input)

        for arg in self.argobs:
            if not arg.display:
                if arg.default is not None:
                    if arg.long_form.startswith('-'):
                        optional_parts.append(arg.long_form)
                        optional_parts.append(str(arg.default))
                    else: # postional args have no name
                        positional_parts.append(str(arg.default))

        # join space separated components of strings with quotes
        for i, part in enumerate(optional_parts):
            if type(part) == str  and ' ' in part:
                optional_parts[i] = "'" + part + "'"

        cmd = ' '.join(optional_parts) + ' ' + ' '.join(positional_parts)
        cmd = cmd.strip()
        print cmd
        return cmd

    def createMouseOverHelp(self, arg):
        """ glue together the arg fields as string """
        help_strings = []

        if arg.required:
            help_strings.append('Required')

        if arg.arg_type is not None:
            help_strings.append('type:' + \
                    str(arg.arg_type).lstrip('<type:class').rstrip('>'))
        elif arg.choices is not None:
            help_strings.append('choice:' + ','.join(map(str, arg.choices)))
        elif arg.action is not None:
            help_strings.append('action:' + str(arg.action))
        else:
            help_strings.append('type: str')

        if arg.default:
            help_strings.append('default: ' + str(arg.default))

        if arg.nargs is not None:
            help_strings.append('nargs: ' + str(arg.nargs))

        if arg.help is not None:
            help_strings.append(arg.help)

        return ', '.join(help_strings)

    # Button press functions

    def onOk(self):
        """ create command-line from form and return it """
        self.accept()

    def onLoad(self):
        """ Create a dialog to open a saved option configuration """
        # TODO: We need a parser to read the config file as a series of ArgObs
        #filename = QtGui.QFileDialog.getOpenFileName()
        #if filename:
        #    helper = argparse.ArgumentParser(add_help=False, parents=[self.parser], fromfile_prefix_chars='@')
        #    self.resetAllWidgets(helper)
        #    result = helper.parse_args(['@{0}'.format(filename)])
        #    for a in helper._get_optional_actions():
        #        self.copyActionValuesToUi(a,result)
        #    for a in helper._get_positional_actions():
        #        self.copyActionValuesToUi(a,result)
        #    self.filename = filename
        pass

    def onSave(self):
        """ Save current option configuration to file """
        # TODO: This is from argparseUI but we actually need to save ArgObs
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
                        "Couldn't write to: {0}".format\
                        (os.path.abspath(self.config_fn)) )

    def onSaveAs(self):
        """ Create a dialog to save option configuration """
        config_fn = QtGui.QFileDialog.getSaveFileName()
        if config_fn:
            QtGui.QMessageBox.information(self,'Information',
                    'Saving config to: {0}'.format\
                    (config_fn))
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
        ArgOb('--group_location', default='all', important=True,
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
        ArgOb('--save_file', type=SaveFilePath, help='path save file'),
        ArgOb('--text', default='rahrahrah', help='whatever')
    ]

    a = AutoGUI(args)
    a.show()
    app.exec_()
    #print ("Ok" if a.result() == 1 else "Cancel")
    if a.result() == 1: # Ok pressed
        print a.makeCommandLine()
        return a.parse_args(a.makeCommandLine())
    else:
        print 'Execution cancelled.'
        sys.exit(0)

if __name__ == '__main__':
    main()