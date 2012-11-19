import argparse
from chippy.ref.util import chroms
from chippy.express import db_query

__author__ = 'Gavin Huttley, Cameron Jack'
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
    def _make_sample_choices(self, db_path):
        """returns the available choices for samples in the DB"""
        session = db_query.make_session('sqlite:///%s' % db_path)
        samples = ['%s : %s' % (s.name, s.description)
               for s in db_query.get_target_sample(session)]
        # These are valid null samples
        samples.insert(0, 'None')
        samples.insert(0, 'none')
        session.close()
        return samples

    def _inc_arg(self, *args, **kwargs):
        """ checks if each potential argument matches
        the user provided list """
        for arg in args:
            tmp_arg = arg.lstrip('-')
            if tmp_arg in self.arg_names:
                if self.required:
                    self.parser.add_argument(*args, required=True, **kwargs)
                else:
                    self.parser.add_argument(*args, **kwargs)

    def add_args(self, arg_names, db_path=None, required=False):
        """ calls _inc_args on every possible argument """

        # Sets whether the options below will be required inputs
        self.required = required
        self.arg_names = arg_names

        self._inc_arg('-c', '--collection', help='Path to the plottable data')
        self._inc_arg('-m', '--metric',
                choices=['Mean counts', 'Frequency counts', 'Standard deviation'],
                default='Frequency counts',
                help='Select the metric (note you will need to change your ylim '\
                'accordingly if providing via --ylim)')
        # Output filename
        self._inc_arg('--plot_filename',
                help='Name of final plot file')
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

        # Sampling choices

        # chrom choice
        self._inc_arg('-C', '--chrom', default='All',
                help='Choose a chromosome',
                choices=('All',)+chroms['mouse'])

        # or external sample (gene) choice
        # NOTE: Should be target_sample
        if db_path:
            self._inc_arg('-E', '--external_sample', default=None,
                    choices=self._make_sample_choices(db_path),
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

        # Misc additional optionals that don't fit an obvious category go here
        self._inc_arg('-p', '--plot_series', action='store_true',
                default=False, help='Plot series of figures. A directory called '\
                +'plot_filename-series will be created.')

        # Not sure what this even does
        self._inc_arg('--text_coords', default=None,
                help='x, y coordinates of series text (e.g. 600,3.0)')

        self._inc_arg('-t', '--test_run', action='store_true',
                help="Test run, don't write output",
                default=False)


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

        self._inc_arg('-s','--sample', default=None,
            choices=self._make_sample_choices(db_path),
            help='Select an existing sample to use')
        self._inc_arg('-S','--new_sample', default=None,
            help="Replace the text on the left and right of the ', "\
            "e.g. `S : S phase'")

        exp_absolute = 'Expression data: absolute ranked'
        exp_diff = 'Expression data: difference in expression between samples'
        target_genes = 'Target gene list'
        self._inc_arg('--sample_type',
                choices=[exp_absolute, exp_diff, target_genes],
                help='Select the type of data you want entered from %s' %\
                str([exp_absolute, exp_diff, target_genes]))

        self._inc_arg('--reffile1', default=None, help='Related file 1')
        self._inc_arg('--reffile2', default=None, help='Related file 2')

    def __init__(self):
        self.parser = argparse.ArgumentParser(description='All ChipPy options')

    # Usage
    # blah = Args()
    # blah.add_args(['collection', 'plot_filename'], required=True)
    # blah.add_args(['ylim', 'xlim'])
    # blah.parser.parse_args()
