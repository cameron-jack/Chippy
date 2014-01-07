"""
    ArgObs is an interface layers between command_args and the
    GUI builder and argument parser of your choosing. Currently
    Argparse and Optparse/CogentOpts are supported.

    Its purpose is to allow more flexibility in passing information
    to the GUI builder, so that more specific interface types can be
    built than is possible with the existing Optparse/OptparseGui
    and Argparse/ArgparseUi module combinations.
"""

__author__ = 'Cameron Jack'
__copyright__ = 'Copyright 2014-2014, Cameron Jack'
__credits__ = ['Cameron Jack']
__license__ = 'GPL'
__version__ = '0.9.dev'
__maintainer__ = 'Cameron Jack'
__email__ = 'cameron.jack@anu.edu.au'
__status__ = 'Pre-release'

import sys
import argparse

from cogent.util.option_parsing import make_option

try:
    from PyQt4 import QtGui
    from gui import Parse2GUI
    GUI_CAPABLE = True
except ImportError:
    print 'Install PyQt4 and ArgparseUi modules to enable GUI'
    GUI_CAPABLE = False

class FilePath(str):
    """
        Holds a string containing the path to a file. Is identifiable as a
        separate type to 'string'. This information is used by the GUI
        builder to create File Browser dialogs and such.
    """
    def __init__(self, path=''):
        super(FilePath, self).__init__()
        self.path = path

class DirPath(str):
    """
        As per file path but dialog will select directory only.
        Use: dir = str(QFileDialog.getExistingDirectory\
                    (self, "Select Directory"))
    """
    def __init__(self, path=''):
        super(DirPath, self).__init__()
        self.path = path

class ArgOb(object):
    """
        An ArgOb is a precursor container for holding command-line
        argument parser entries. This exists for two reasons:
        1) You can't iterate through Argparse parser arguments
        2) Extra information than what is typically allowed in
        Argparse argument construction can be provided to automatic
        GUI builders that in-turn create the command-line that gets
        parsed.

        An example of 2. is the custom type FilePath provided here,
        another is the 'display' field, which allows GUI building to
        be turned off for an arg, though it will still be created.
        This allows a 'pass-thru' capability for predefined options
        like which DB to use.
    """
    def __init__(self, *args, **kwargs):
        """
            *args are short-form and long-form
            **kwargs are nargs, type, action, help, default, display
            'display' defaults to True but if False will be passed to the
            command line but not displayed in the GUI.
        """
        if len(args) == 1:
            if args[0].startswith('-'):
                self.positional = False
            else:
                self.positional = True
            self.short_form = None
            self.long_form = args[0]
        elif len(args) == 2:
            if args[0].startswith('--'):
                # this shouldn't happen, but just in case...
                self.short_form = args[1]
                self.long_form = args[0]
            else:
                self.short_form = args[0]
                self.long_form = args[1]
            self.positional = False

        self.arg_type = None
        if 'type' in kwargs.keys():
            self.arg_type = kwargs['type']
        self.nargs = None
        if 'nargs' in kwargs.keys():
            self.nargs = kwargs['nargs']
        self.action = None
        if 'action' in kwargs.keys():
            self.action = kwargs['action']
        self.choices = None
        if 'choices' in kwargs.keys():
            self.choices = kwargs['choices']

        self.required = False
        if 'required' in kwargs.keys():
            self.required = kwargs['required']

        self.default = None
        if 'default' in kwargs.keys():
            self.default = kwargs['default']
        self.metavar = None
        if 'metavar' in kwargs.keys():
            self.metavar = kwargs['metavar']

        self.display = True
        if 'display' in kwargs.keys():
            self.display = kwargs['display']

        self.help = ''
        if 'help' in kwargs.keys():
            self.help = kwargs['help']

        self.cogent_type = None
        if 'cogent_type' in kwargs.keys():
            self.cogent_type = kwargs['cogent_type']

    _ARG_TO_OPT_TYPE_CONVERTER = {None: 'string', int: 'int',
            long: 'long', float: 'float', 'choices': 'choice',
            FilePath: 'existing_filepath', DirPath: 'existing_dirpath'}

    def asCogentOpt(self):
        """
            Creates and returns a PyCogent CogentOption object from
            an argparse entry. Called by _inc_arg().
            Also creates entries for the optparse_gui graphical interface.

            The following types need to available for full type conversion to
            Galaxy type to take place. Types 'string' through 'choice' are
            implied by their argparse type entry and do NOT need to be
            specified explicitly. The types from 'multiple_choice' through
            'new_path' DO need to be explicitly provided as 'cogent_type'
            entries.

            CogentOption.TYPE on left, Galaxy_type on right:

            type_converter['string'] = "text"
            type_converter['int'] = "integer"
            type_converter['long'] = "float"
            type_converter['float'] = "float"
            type_converter['choice'] = "select"

            type_converter['multiple_choice'] = "multiple_select"
            type_converter['existing_filepath'] = "data"
            type_converter['existing_filepaths'] = "repeat"
            type_converter['existing_dirpath'] = "input_dir"
            type_converter['existing_path'] = "input_dir"
            type_converter['new_filepath'] = "output"
            type_converter['new_dirpath'] = "output_dir"
            type_converter['new_path'] = "output_dir"
        """

        cogent_parameters = {}
        if self.action is not None:
            cogent_parameters['action'] = self.action
        elif self.cogent_type is not None:
            cogent_parameters['cogent_type'] = self.cogent_type
        elif self.choices is not None:
            cogent_parameters['choices'] = self.choices
            cogent_parameters['type'] = 'choice'
        elif self.arg_type is not None:
            try:
                cogent_parameters['type'] =\
                        self._ARG_TO_OPT_TYPE_CONVERTER[self.arg_type]
            except KeyError:
                cogent_parameters['type'] = 'string'
        cogent_parameters['help'] = self.help

        if self.short_form is not None:
            return make_option(self.short_form, self.long_form,
                    cogent_parameters)
        else:
            return make_option(self.long_form, **cogent_parameters)

    def addToArgparse(self, parser):
        """ create an argparse arg and add it to the parser """
        argparse_params = {}

        # nargs, metavar, help and required are optional
        if self.nargs is not None:
            argparse_params['nargs'] = self.nargs
        if self.metavar is not None:
            argparse_params['metavar'] = self.metavar
        if self.help is not None:
            argparse_params['help'] = self.help
        if self.required is not None:
            argparse_params['required'] = self.required

        # type, choices and action are mutually exclusive
        if self.arg_type is not None:
            argparse_params['type'] = self.arg_type
        elif self.choices is not None:
            argparse_params['choices'] = self.choices
        elif self.action is not None:
            argparse_params['action'] = self.action

        if self.short_form is not None:
            parser.add_argument(self.short_form, self.short_form,
                    argparse_params)
        else:
            parser.add_argument(self.long_form, **argparse_params)

def test_argobs():
    """ Create some argobs and add them to a parser """
    parser = argparse.ArgumentParser()

    arg1 = ArgOb('--num_genes', type=int, required=True,
        help='Number of ranked genes to use in study')

    arg2 = ArgOb('--group_location', default='all',
        choices=['all', 'top', 'middle', 'bottom'],
        help='The representative group in a study to form a plot line')

    arg3 = ArgOb('--ranks', action='store_true',
        help='Use rank-based data values instead of counts or '+\
             'expression')

    arg4 = ArgOb('--sample_extremes', type=float, default=0.0,
            help='Proportion of least and most absolute '+\
            'expressed genes to treat separately. Set to 0.0 to disable.')

    arg1.addToArgparse(parser)
    arg2.addToArgparse(parser)
    arg3.addToArgparse(parser)
    arg4.addToArgparse(parser)

    cmd_line = '--num_genes 4 --group_location all --ranks --sample_extremes 0.4'.split(' ')
    args = parser.parse_args(cmd_line)
    print args.__dict__

if __name__ == '__main__':
    test_argobs()


