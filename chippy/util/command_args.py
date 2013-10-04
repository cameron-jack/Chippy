import argparse
import sys # required for finding db_path
from chippy.express.db_query import get_chroms
from chippy.express import db_query
from chippy.express.util import sample_types

from cogent.util.option_parsing import make_option
try:
    from PyQt4 import QtGui
    from argparseui.argparseui import ArgparseUi
    GUI_CAPABLE = True
except ImportError:
    print 'Install PyQt4 and ArgparseUi modules to enable GUI'
    GUI_CAPABLE = False

__author__ = 'Cameron Jack'
__copyright__ = 'Copyright 2011-2013, Gavin Huttley, Anuj Pahwa, Cameron Jack'
__credits__ = ['Cameron Jack']
__license__ = 'GPL'
__maintainer__ = 'Cameron Jack'
__email__ = 'cameron.jack@anu.edu.au'
__status__ = 'pre-release'
__version__ = '638'

"""
    command_args offers the entire arguments/options set to define the
    interfaces for all ChipPy scripts. It returns both an argparse parser
    object and a script_info{} entry with PyCogent CogentOption objects
    so that the Qiime xml_generator can be used to automagically create
    Galaxy interfaces for each script.

    The following types need to available for full type conversion to
    Galaxy type to take place. Types 'string' through 'choice' are implied
    by their argparse type entry and do NOT need to be specified explicitly.
    The types from 'multiple_choice' through 'new_path' DO need to be
    explicitly provided as 'cogent_type' entries.

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

# Example of mutually exclusive options
# group = parser.add_mutually_exclusive_group()
# Use required=True to require an argument

class Args(object):
    def parse(self):
        """
            Returns the result of calling the parser on any input.
            Supports both command line and graphic interface through
            ArgparseUi and PyQt4.
        """
        # This will either use argparse if arguments are provided or will
        # call the ArgparseUi graphic interface (PyQT)
        if GUI_CAPABLE:
            app = QtGui.QApplication(sys.argv)
            a = ArgparseUi(self.parser)
            a.show()
            app.exec_()
            print ("Ok" if a.result() == 1 else "Cancel")
            if a.result() == 1: # Ok pressed
                return a.parse_args()
            else:
                return self.parser.parse_args()
        else:
            return self.parser.parse_args()

    def _make_sample_choices(self):
        """returns the available choices for samples in the DB"""
        if self.db_path is None:
            return ['None', 'none']
        session = db_query.make_session(str(self.db_path))
        samples = db_query.get_sample_choices(session)
        #samples = ['%s : %s' % (str(s.name), str(s.description))
        #        for s in db_query.get_samples(session)]
        # These are valid null samples
        #samples.insert(-1, 'None')
        #samples.insert(-1, 'none')
        session.close()
        return samples

    def _make_chrom_choices(self):
        if self.db_path is None:
            return []
        session = db_query.make_session(str(self.db_path))
        chroms = get_chroms(session)
        if not chroms: # must return something or opt build will fail
            chroms = ['None', 'none']
        return chroms

    # argparse base arg type converter to optparse basic types allows us to
    # avoid specifically adding 'cogent_action_or_type' to every argument.
    _ARG_TO_OPT_TYPE_CONVERTER = {None: "string", int: "int",
            long: "long", float: "float", 'choices': 'choice'}

    def _make_cogent_opt(self, cogent_type, *args, **kwargs):
        """
            Creates and returns a PyCogent CogentOption object from
            an argparse entry. Called by _inc_arg().
            Also creates entries for the optparse_gui graphical interface
        """
        if 'choices' in kwargs:
            kwargs['type'] = 'choice'
        elif 'action' in kwargs:
            # Required, because without this, type=None implies type='string'
            pass
        else:
            if cogent_type:
                kwargs['type'] = cogent_type
            else:
                arg_type = kwargs.get('type', None)
                kwargs['type'] = self._ARG_TO_OPT_TYPE_CONVERTER[arg_type]
        if len(args) == 2:
            cogent_opt = make_option(args[0], args[1], **kwargs)
        else:
            cogent_opt = make_option(args[0], **kwargs)
        return cogent_opt

    def _inc_arg(self, *args, **kwargs):
        """ Checks if each potential argument matches the user provided list.
            Then creates an argparse parser entry and a PyCogent
            CogentOption object. """

        # pull out cogent_type now so that add_argument still works
        cogent_type = kwargs.pop('cogent_type', None)
        for arg in args:
            tmp_arg = arg.lstrip('-')
            # Required args
            if self.required_args and tmp_arg in self.required_args:
                # add argparse entry
                self.parser.add_argument(*args, required=True, **kwargs)
                # add CogentOption entry
                cogent_opt = self._make_cogent_opt(cogent_type, *args, **kwargs)
                self.req_cogent_opts.append(cogent_opt)

            # Optional args
            elif self.optional_args and tmp_arg in self.optional_args:
                # add argparse entry
                self.parser.add_argument(*args, **kwargs)
                # add CogentOption entry
                cogent_opt = self._make_cogent_opt(cogent_type, *args, **kwargs)
                self.opt_cogent_opts.append(cogent_opt)


    def _add_load_save_args(self):
        """ All loading and saving related arguments should go here """

        self._inc_arg('--make_bedgraph', action='store_true', help='Enable '+\
                'Output to BEDgraph during export. Save name is same as ChipPy DB'+\
                "but with a '_expression-name.bedgraph' extension")
        self._inc_arg('-c', '--collection', help='Path to the plottable data')
        self._inc_arg('--plot_filename',
                help='Name of final plot file')
        self._inc_arg('-e', '--expression_data',
                help="Path to the expression/gene data file. Must be tab delimited.")
        self._inc_arg('--allow_probeset_many_gene', action='store_true',
                default=False, help='Allow probesets that map to multiple genes')
        self._inc_arg('--gene_id_heading',
                default='gene',
                help='Column containing the Ensembl gene stable ID')
        self._inc_arg('--probeset_heading',
                default='probeset',
                help='Column containing the probeset IDs')
        self._inc_arg('--expression_heading',
                default='exp',
                help='Column containing the expression scores')

        # Only valid if a CHIPPY_DB has been given
        if self.db_path:
            self._inc_arg('-s','--sample', choices=self._make_sample_choices(),
                    help="Select an existing sample to use, form is '+\
                    'S : S phase'")
            # Generic multi-sample options (perhaps should be positional)
            self._inc_arg('--sample1', choices=self._make_sample_choices(),
                    help="Choose a first expression study, form is '+\
                    'S : S phase'")
            self._inc_arg('--sample2', choices=self._make_sample_choices(),
                    help="Choose a second expression study, form is '+\
                    'S : S phase'")
            self._inc_arg('--sample3', choices=self._make_sample_choices(),
                    help="Choose a third expression study, form is '+\
                    'S : S phase'")
            self._inc_arg('--sample4', choices=self._make_sample_choices(),
                    help="Choose a fourth expression study, form is '+\
                    'S : S phase'")
            samplesStr = ', '.join(self._make_sample_choices())
            self._inc_arg('-S','--new_sample', default=None,
                    help="Select an existing sample to use, form is '+\
                    'S : S phase'. Existing samples: "+samplesStr)
            # It would be nice if this only offered diff studies as choices
            self._inc_arg('--diff_sample', help='Choose the expression study',
                    choices=[str(s) for s in self._make_sample_choices()])

        self._inc_arg('--sample_type',
                choices=[k for k in sample_types.keys()],
                help='Select the type of data you want entered from '+\
                ' '.join([str(k) for k in sample_types.keys()]) )

        self._inc_arg('--reffile1', default=None, help='Related file 1')
        self._inc_arg('--reffile2', default=None, help='Related file 2')

        # Export Centred Counts args
        self._inc_arg('-B', '--BAMorBED',
                help='Read counts. Either an indexed BAM or a BED file')

        self._inc_arg('-f', '--overwrite', action='store_true',
                help='Ignore any saved files', default=False)
        self._inc_arg('--tab_delimited', action='store_true',
                help='output to tab delimited format', default=False)
        self._inc_arg('--chr_prefix', default='', help='String added by '+\
                'aligner to prefix chromosome numbers/names.')

    def _add_sampling_args(self):
        """ All arguments relate to conditional selection of data """
        # chrom choice
        self._inc_arg('-C', '--chrom', default=None,
                help='Choose a chromosome',
                choices=(self._make_chrom_choices()) )

        # group genes into sets ranked by expression
        self._inc_arg('-g', '--group_size', default='All',
                help='Number of genes to group to estimate'+\
                 ' statistic - All or a specific number')

        self._inc_arg('--num_genes', type=int,
                help='Number of ranked genes to use in study')

        self._inc_arg('--group_location', default='all',
                choices=['all', 'top', 'middle', 'bottom'],
                help='The representative group in a study to form a plot line')

        self._inc_arg('--ranks', action='store_true',
                help='Use rank-based data values instead of counts or '+\
                'expression')

        self._inc_arg('--sample_extremes', type=float,
                default=0.0, help='Proportion of least and most absolute '+\
                'expressed genes to treat separately. Set to 0.0 to disable.')

        # plot only the most expressed genes of group_size
        # DEPRECATED
        self._inc_arg('--top_features', action='store_true', default = False,
                help='Plot only top features ranked by expressed chromatin')

        # Use moving-mean (boxcar) smoothing on lines
        self._inc_arg('--smoothing', type=int, default=0,
                help='Window size for smoothing of plot data')

        self._inc_arg('--binning', type=int, default=0,
                help='Sum counts within integer-sized bins across plot')

        # Filter out over- and under-expression outliers
        self._inc_arg('-k', '--cutoff', type=float, default = 0.05,
                help='Probability cutoff. Exclude genes if the probability of '+\
                'the observed tag count is less than or equal to this '+\
                'value e.g. 0.05.')

        # Export Centred Counts args
        self._inc_arg('--expression_area', choices=['TSS', 'Exon_3p',
                'Intron_3p', 'Both_3p'], help='Expression area options: '+\
                'TSS, Exon_3p, Intron-3p, Both-3p')

        self._inc_arg('--window_radius', type=int, default=1000,
                help='Region size around TSS')

        self._inc_arg('--multitest_signif_val', type=int,
                help='Restrict plot to genes that pass multi-test '+\
                'significance. Valid values: 1, 0, -1', default=None)

        self._inc_arg('--include_target', default=None,
                help='A Target Gene List in ChipPyDB',
                choices=[str(s) for s in self._make_sample_choices()])
        self._inc_arg('--exclude_target', default=None,
                help='A Target Gene List in ChipPyDB',
                choices=[str(s) for s in self._make_sample_choices()])

    def _add_plot_args(self):
        """ Arguments specifically related to showing graphical plots """
        self._inc_arg('--plot_format', default='png', choices=['png', 'pdf'],
                help="Select the plot format to output: 'png' or 'pdf'"+\
                "[default: %default]")

        self._inc_arg('-y', '--ylim', default=None,
                help='comma separated minimum-maximum yaxis values (e.g. 0,3.5)')
        self._inc_arg('-H', '--fig_height', type=float, default=6*2.5,
                help='Figure height (cm)')
        self._inc_arg('-W', '--fig_width', type=float, default=10*2.5,
                help='Figure width (cm)')

        # Important note, grid_lines are an absolute scale!
        self._inc_arg('--xgrid_lines', type=float, default = 100,
                help='major grid-line spacing on x-axis')
        self._inc_arg('--ygrid_lines', type=float, default = None,
                help='major grid-line spacing on y-axis')
        self._inc_arg('--grid_off', action='store_true',
                help='Turn grid lines off')
        # tick spacing along axes
        self._inc_arg('--xtick_interval', type=int, default=2,
                help='number of blank ticks between labels')
        self._inc_arg('--ytick_interval', type=int, default=2,
                help='number of blank ticks between labels')
        # Smooth top and right borders for cleaner looking plot
        self._inc_arg('--clean_plot', action='store_true', default=False,
                help='Remove tick marks and top and right borders ')
        # background colour - black or white
        self._inc_arg('-b', '--bgcolor', default='black',
                help='Plot background color',
                choices=['black', 'white'])
        # side colour bar for 3rd-dimension range
        self._inc_arg('--colorbar', action='store_true',
                help="Add colorbar to figure", default=False)

        # Create options for font sizes, colours and labels
        # Headline of the plot
        self._inc_arg('--title', help='Plot title')
        # Need options for size and style
        # Axis labels and font sizes
        self._inc_arg('--ylabel', default = 'Normalized counts',
             help='Label for the y-axis')
        self._inc_arg('--xlabel', default = 'Position relative to TSS',
                help='Label for the x-axis')
        self._inc_arg('--xfont_size', type=int, default=12,
                help='font size for x label')
        self._inc_arg('--yfont_size', type=int, default=12,
                help='font size for y label')
        # Optional axis units text
        self._inc_arg('--yaxis_units',
                help='Text showing units of y-axis of plot')
        self._inc_arg('--xaxis_units',
                help='Text showing units of x-axis of plot')

        # Turn on line legend and set font size
        self._inc_arg('-l', '--legend', action='store_true', default=False,
                help='Automatically generate a figure legend.')
        self._inc_arg('--legend_font_size', type=int, default=12,
                help='Point size for legend characters')

        # Vertical line showing the centred feature
        self._inc_arg('--vline_style',
                default = '-.', choices=['-.', '-', '.'],
                help='line style for centred vertical line')
        self._inc_arg('--vline_width', type=int,
                default = 2,
                help='line width for centred vertical line')

        # Plotted line opacity
        self._inc_arg('--line_alpha', type=float, default=1.0,
                help='Opacity of lines')

        # Plot as grey-scale
        self._inc_arg('--grey_scale', action='store_true',
                help='Plot colour range as grey-scale')

        # Plot as series of images
        self._inc_arg('--plot_series', action='store_true',
                help='Plot series of figures. A directory called '+\
                'plot_filename-series will be created.')

        # Position of legend?
        self._inc_arg('--text_coords', default=None,
                help='x, y coordinates of series text (e.g. 600,3.0)')

        self._inc_arg('--div', help='Path to the plottable data, to divide '+\
                'other loaded plottable data')

        self._inc_arg('--confidence_intervals', action='store_true',
                help='Show confidence intervals around plot lines')

        # Uses the 'tag count' arg value in the export file
        self._inc_arg('--normalise_by_RPM', action='store_true',
                help='Normalise by Reads Per mapped-Million')

    def _add_db_args(self):
        """ options for starting a ChipPy DB """
        self._inc_arg('--save_db_path',
                help='path to directory where chippy.db will be saved')
        self._inc_arg('--ensembl_release', type=int,
                help='Ensembl release to use.')
        self._inc_arg('--species', default='mouse',
                help="Create for species e.g. 'mouse','human'")
        self._inc_arg('--hostname', default=None,
                help='hostname for SQL Ensembl server')
        self._inc_arg('--username', default=None,
                help='username SQL Ensembl server')
        self._inc_arg('--password', default=None,
                help='password for SQL Ensembl server')
        self._inc_arg('--port', default=None, type=int,
                help='Port for SQL Ensembl server')
        self._inc_arg('--dummy_data', default=False, action='store_true',
                help='Create dummy expression data sets')

    def _add_misc_args(self):
        """ various options that don't fall into a category above """

        self._inc_arg('-m', '--metric',
                choices=['Mean counts', 'Frequency counts', 'Standard deviation'],
                default='Frequency counts',
                help='Select the metric (note you will need to change your ylim '\
                'accordingly if providing via --ylim)')

        self._inc_arg('-t', '--test_run', action='store_true',
                help="Test run, don't write output",
                default=False)

    def _add_diff_abs_plots_specific_args(self):
        """ These args are specifically for the diff_abs_plots script """

        # dot colouring
        self._inc_arg('--extremes_colour', default='blue', choices=['blue',
                'red', 'yellow', 'green', 'magenta', 'orange', 'cyan'],
                help='Colour of dots for absolute expression marked '\
                'as extreme.')
        self._inc_arg('--signif_colour', default='blue', choices=['blue',
                'red', 'yellow', 'green', 'magenta', 'orange', 'cyan'],
                help='Colour of dots for difference of expression marked '\
                'as significant.')
        self._inc_arg('--bulk_colour', default='blue', choices=['blue',
                'red', 'yellow', 'green', 'magenta', 'orange', 'cyan'],
                help='Colour of dots for all relatively unexceptional '\
                'expression values.')

        # hide unwanted plot areas
        self._inc_arg('--hide_extremes', default=False, action='store_true',
                help='Do not show absolute expression considered extreme')
        self._inc_arg('--hide_signif', default=False, action='store_true',
                help='Do not show difference expression considered '+\
                'significant')
        self._inc_arg('--hide_bulk', default=False, action='store_true',
                help='Do not show expression values considered normal')

        # misc. options that don't fit other categories
        self._inc_arg('--gene_file', default=None,
                help='Annotated gene list file output path, as pickle.gz')
        self._inc_arg('--output_prefix1',
                default=None, help='Output path prefix for first plot')
        self._inc_arg('--output_prefix2',
                default=None, help='Output path prefix for second plot')

    def _add_ChrmVsExpr_specific_args(self):
        """ These args are specifically for the ChrmVsExpr script """
        self._inc_arg('--dot_or_line_plot', choices=['dot', 'line'],
                help='Select the type of plot you want '\
                'entered from %choices')
        self._inc_arg('--x_axis_type',
                choices=['expression', 'chrm counts', 'expr counts'],
                help='Select the data unit type for the x-axis')
        self._inc_arg('--y_axis_type',
                choices=['expression', 'chrm counts', 'expr counts'],
                help='Select the data unit type for the y-axis')


    def __init__(self, positional_args=None, required_args=None,
            optional_args=None):
        """ calls _inc_args on every possible argument """
        self.parser = argparse.ArgumentParser(description=\
                'All ChipPy options', version='ChipPy r'+str(__version__))
        # optgui_parser is the optparse_gui parser object
        self.req_cogent_opts = []
        self.opt_cogent_opts = []
        self.db_path = None

        # We need to handle the 'db_path' positional argument ourselves
        # as several arguments require it already loaded
        if positional_args:
            if 'db_path' in positional_args:
            # Need to add an actual positional arg now!
                self.parser.add_argument('db_path',
                        help='Path to ChippyDB')
                cogent_opt = make_option('--db_path',
                        type='existing_filepath', help='Path to ChippyDB')
                self.req_cogent_opts.append(cogent_opt)
                for arg in sys.argv:
                    if not arg.startswith('-'):
                        # it's positional so take this to be db_path
                        possible_db_path = str(arg).strip()
                        if 'chippy' in possible_db_path.lower() and\
                                '.db' in possible_db_path.lower():
                            self.db_path = possible_db_path
                            print self.db_path, 'selected as ChipPy database'
                            break

        self.required_args = required_args
        self.optional_args = optional_args

        # process arguments for loading or saving
        self._add_load_save_args()
        # process arguments for sub-sampling
        self._add_sampling_args()
        # process arguments for graphical plotting
        self._add_plot_args()
        # process arguments for creating the DB
        self._add_db_args()
        # process misc arguments
        self._add_misc_args()
        # process args specific to the diff_abs_plots script
        self._add_diff_abs_plots_specific_args()
        # process args specific to the ChrmVsExpr script

