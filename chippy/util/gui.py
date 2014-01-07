import textwrap
import argparse
from PyQt4 import QtCore, QtGui

from command_args import FilePath

__author__ = 'Cameron Jack'
__copyright__ = 'Copyright 2014-2014, Cameron Jack'
__credits__ = ['Cameron Jack']
__license__ = 'GPL'
__version__ = '0.9.dev'
__maintainer__ = 'Cameron Jack'
__email__ = 'cameron.jack@anu.edu.au'
__status__ = 'Pre-release'

class HelpText(object):
    """ Construction of argument information to text help for display """
    def __init__(self, long_form, num_args, short_form = '', type_or_action='', is_required='Not required'):
        self.num_args = num_args
        self.type_or_action = type_or_action
        self.is_required = is_required
        self.short_form = short_form
        self.long_form = long_form

    def constructHelpText(self):
        """ based on available info, join fields to form sentences """
        parts = [self.short_form,
                 self.long_form,
                 self.num_args,
                 self.type_or_action,
                 self.is_required
        ]
        return ' '.join(parts)

class Parse2GUI(QtGui.QDialog):
    """
        Interprets Argparse arguments and builds a GUI for the user,
        that is then parsed by Argparse and returned.

        Specifically uses GridLayout to allow more flexibility with
        widget types and interface design - in particular the use of
        buttons to launch dialog boxes.
    """

    def __init__(self, parser, window_title='Make your choice', parent=None):
        super(Parse2GUI, self).__init__(parent)
        self.setWindowTitle(window_title)
        self.parser = parser

        self.createOptions()

    def createOptions(self):
        """
            Create option names (long form only) and their input field
            as appropriate for their type and action. Files should open
            a browser dialog so that the user can easily select the
            locations.
        """

        # fp = FilePath()
        # for i, arg in enumerate(arguments):
        #     if arg.type == type(fp):
        #         # Create file dialog button and lineEdit widget pair
        #     elif arg.type == int:
        #         try:
        #             nargs = int(args.nargs)
        #         except ValueError:
        #             nargs = -1
        #         if arg.nargs == '*':
        #
        name_line_edit = QtGui.QLineEdit()
        name_label = QtGui.QLabel('File1:')
        name_label.setBuddy(name_line_edit)

        email_line_edit = QtGui.QLineEdit()
        email_label = QtGui.QLabel('Name:')
        email_label.setBuddy(email_line_edit)

        age_spin_box = QtGui.QSpinBox()
        age_label = QtGui.QLabel('Name:')
        age_label.setBuddy(age_spin_box)

        gridLayout = QtGui.QGridLayout()
        gridLayout.addWidget(name_label, 0, 0)
        gridLayout.addWidget(name_line_edit, 0, 1)
        gridLayout.addWidget(email_label, 1, 0)
        gridLayout.addWidget(email_line_edit, 1, 1)
        gridLayout.addWidget(age_label, 2, 0)
        gridLayout.addWidget(age_spin_box, 2, 1)
        QtGui.setLayout(gridLayout)

        self.OkButton = self.addButton("Ok")
        self.CancelButton= self.addButton("Cancel")
        self.OkButton.clicked.connect(self.onOk)
        self.CancelButton.clicked.connect(self.onCancel)

    def addButton(self, label):
        self.buttonsLayout.addSpacerItem(QtGui.QSpacerItem(\
                20, 1, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
        b = QtGui.QPushButton(label, self.buttons)
        self.buttonsLayout.addWidget(b)
        return b

    def onOk(self):
        """ Ok button pressed """
        self.accept()

    def onCancel(self):
        """ Cancel button pressed """
        self.reject()