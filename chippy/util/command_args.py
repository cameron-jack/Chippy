import argparse
import sys # required for finding db_path
from chippy.ref.util import chroms
from chippy.express import db_query
from chippy.express.util import sample_types

__author__ = 'Cameron Jack'
__copyright__ = 'Copyright 2011, Gavin Huttley, Anuj Pahwa, Cameron Jack'
__credits__ = ['Gavin Huttley, Cameron Jack']
__license__ = 'GPL'
__maintainer__ = 'Gavin Huttley'
__email__ = 'Gavin.Huttley@anu.edu.au'
__status__ = 'alpha'
__version__ = '0.1'

# Example of mutually exclusive options
#group = parser.add_mutually_exclusive_group()
# Use nargs='1' or nargs='+' to require 1 or, 1 or more arguments
# Use required=True to require an argument

class Args(object):
    def _make_sample_choices(self):
        """returns the available choices for samples in the DB"""
        if self.db_path is None:
            return ['None', 'none']
        session = db_query.make_session('sqlite:///' + str(self.db_path))
        samples = ['%s : %s' % (str(s.name), str(s.description))
                for s in db_query.get_samples(session)]
        # These are valid null samples
        samples.insert(-1, 'None')
        samples.insert(-1, 'none')
        session.close()
        return samples

    def _inc_arg(self, *args, **kwargs):
        """ checks if each potential argument matches
        the user provided list """
        for arg in args:
            tmp_arg = arg.lstrip('-')
            if self.required_args and tmp_arg in self.required_args:
                self.parser.add_argument(*args, required=True, **kwargs)
            elif self.optional_args and tmp_arg in self.optional_args:
                self.parser.add_argument(*args, **kwargs)

    def _add_load_save_args(self):
        """ All loading and saving related arguments should go here """

        # Plot centred counts args
        self._inc_arg('-c', '--collection', help='Path to the plottable data')
        # Output filename
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

        if self.db_path:
            self._inc_arg('-s','--sample', default=None,
                    choices=self._make_sample_choices(),
                    help="Select an existing sample to use, form is '\
                    'S : S phase'")
            samplesStr = ', '.join(self._make_sample_choices())
            self._inc_arg('-S','--new_sample', default=None,
                    help="Select an existing sample to use, form is '\
                    'S : S phase'. Existing choices: "+samplesStr)

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

    def _add_sampling_args(self):
        """ All arguments relate to conditional selection of data """
        # chrom choice
        self._inc_arg('-C', '--chrom', default='All',
            help='Choose a chromosome',
            choices=('All',)+chroms['mouse'])

        # or external sample (gene) choice
        # NOTE: Should be target_sample
        if self.db_path:
            self._inc_arg('-E', '--external_sample', default=None,
                choices=self._make_sample_choices(),
                help='External sample')

        # group genes into sets ranked by expression
        self._inc_arg('-g', '--group_size', default='All',
            help='Number of genes to group to estimate'\
                 ' statistic - All or a specific number')

        self._inc_arg('--group_location', default='all',
            choices=['all', 'top', 'middle', 'bottom'],
            help='The representative group in a study to form a plot line')

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
        self._inc_arg('-k', '--cutoff', type=float, default = None,
            help='Probability cutoff. Exclude genes if the probability of '\
                 'the observed tag count is less than or equal to this value '\
                 'e.g. 0.05.')

        # Export Centred Counts args
        self._inc_arg('--expression_area', choices=['TSS', 'Exon_3p',
                'Intron_3p', 'Both_3p'], help='Expression area options: ' \
                'TSS, Exon_3p, Intron-3p, Both-3p')

        self._inc_arg('--max_read_length', type=int,
            default=75, help='Maximum sequence read length')

        self._inc_arg('--count_max_length',
                action='store_false', help='Use maximum read length instead '\
                          'of mapped length', default=True)

        self._inc_arg('--window_size', type=int, default=1000,
                help='Region size around TSS')

        self._inc_arg('--multitest_signif_val', type=int,
                help='Restrict plot to genes that pass multitest signficance,'\
                'valid values: 1, 0, -1', default=None)

        self._inc_arg('--include_target', default=None,
                help='A Target Gene List in ChipPyDB',
                choices=[str(s) for s in self._make_sample_choices()])
        self._inc_arg('--exclude_target', default=None,
                help='A Target Gene List in ChipPyDB',
                choices=[str(s) for s in self._make_sample_choices()])

    def _add_plot_args(self):
        """ Arguments specifically related to showing graphical plots """
        # Create options essential for making a plot
        self._inc_arg('-y', '--ylim', default=None,
            help='comma separated minimum-maximum yaxis values (e.g. 0,3.5)')
        self._inc_arg('-H', '--fig_height', type=float, default=2.5*3,
            help='Figure height (cm)')
        self._inc_arg('-W', '--fig_width', type=float, default=2.5*5,
            help='Figure width (cm)')

        # Important note, grid_lines are an absolute scale!
        self._inc_arg('--xgrid_lines', type=float, default = 100,
            help='major grid-line spacing on x-axis')
        self._inc_arg('--ygrid_lines', type=float, default = None,
            help='major grid-line spacing on y-axis')
        self._inc_arg('--grid_off', action='store_true', default=False,
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
        self._inc_arg('--x_font_size', type=int, default=12,
            help='font size for x label')
        self._inc_arg('--y_font_size', type=int, default=12,
            help='font size for y label')

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

        self._inc_arg('-p', '--plot_series', action='store_true',
            default=False, help='Plot series of figures. A directory called '\
            +'plot_filename-series will be created.')

        # Position of legend?
        self._inc_arg('--text_coords', default=None,
            help='x, y coordinates of series text (e.g. 600,3.0)')

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

    def __init__(self, positional_args=None, required_args=None,
                optional_args=None):
        """ calls _inc_args on every possible argument """
        self.parser = argparse.ArgumentParser(description=\
                'All ChipPy options')
        self.db_path = None

        # We need to handle the 'db_path' positional argument ourselves
        # as several arguments require it already loaded
        if positional_args:
            if 'db_path' in positional_args:
                for arg in sys.argv:
                    if not '-' in arg:
                        # it's positional so take this to be db_path
                        possible_db_path = str(arg).strip()
                        if 'chippy' in possible_db_path.lower() and\
                                '.db' in possible_db_path.lower():
                            self.db_path = possible_db_path
                            print self.db_path, 'selected as ChipPy database'
                            # Need to add an actual positional arg now!
                            self.parser.add_argument('db_path')
                            break
                if self.db_path is None:
                    raise RuntimeError('Path to ChipPy database required')

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

        self.parsed_args = self.parser.parse_args()
