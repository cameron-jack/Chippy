import argparse
import sys
sys.path.extend(['..','../..'])

from chippy.express import db_query
from chippy.express.util import sample_types, sample_type_names, \
    sample_type_desc

from argobs import OpenFilePath, SaveFilePath, DirPath, ArgOb

try:
    from PyQt4 import QtGui
    from gui import AutoGUI
    GUI_CAPABLE = True
except ImportError:
    print 'Install PyQt4 to enable GUI'
    GUI_CAPABLE = False

__author__ = 'Cameron Jack'
__copyright__ = 'Copyright 2011-2013, Gavin Huttley, Anuj Pahwa, Cameron Jack'
__credits__ = ['Cameron Jack']
__license__ = 'GPL'
__maintainer__ = 'Cameron Jack'
__email__ = 'cameron.jack@anu.edu.au'
__status__ = 'pre-release'
__version__ = '711'

"""
    command_args offers the entire arguments/options set to define the
    interfaces for all ChipPy scripts. It returns both an argparse parser
    object and a script_info{} entry with PyCogent CogentOption objects
    so that the Qiime xml_generator can be used to automagically create
    Galaxy interfaces for each script.
"""

# Example of mutually exclusive options
# group = parser.add_mutually_exclusive_group()
# Use required=True to require an argument

class Args(object):
    def parse(self, **kwargs):
        """
            Returns the result of calling the parser on any input.
            Supports both command line and graphic interface through
            ArgparseUi and PyQt4.
        """
        # This will either use argparse if arguments are provided or will
        # call the ArgparseUi graphic interface (PyQT)
        for argob in self.argobs:
            argob.addToArgparse(self.parser)

        if GUI_CAPABLE:
            if len(sys.argv) < 3:
                app = QtGui.QApplication(sys.argv)
                # Use ArgparseUi or AutoGUI
                #a = ArgparseUi(self.parser, **kwargs)
                a = AutoGUI(self.argobs, **kwargs)
                a.show()
                app.exec_()
                if a.result() == 1: # Ok pressed
                    return self.parser.parse_args(a.makeCommands())
                else:
                    print 'Execution cancelled.'
                    sys.exit(0)
            else:
                return self.parser.parse_args()
        else:
            return self.parser.parse_args()

    def _create_argob(self, *args, **kwargs):
        """
            Checks if each potential argument matches the user provided list.
            Then creates an ArgOb
        """
        for arg in args:
            tmp_arg = arg.lstrip('-')
            # Required args
            if self.required_args and tmp_arg in self.required_args:
                self.argobs.append(ArgOb(*args, required=True, **kwargs))
            elif self.optional_args and tmp_arg in self.optional_args:
                self.argobs.append(ArgOb(*args, **kwargs))

    def _add_load_save_args(self):
        """ All loading and saving related arguments should go here """

        self._create_argob('--make_bedgraph', action='store_true',
                help='Save exported windows to BEDgraph. The output file '+
                'name is the same as the collection file but with an '+\
                "'_expression-name.bedgraph' extension")
        self._create_argob('--BED_windows', action='store_true',
                help='Save region of interest windows to BED. The output '+\
                'file name is the same as collection file but with '+\
                "a '_.BED' extension. Used for annotation.")
        self._create_argob('-c', '--collection', type=SaveFilePath,
                help='Path to save plottable data. These files are given '+\
                "the '.chp' extension for clarity")
        self._create_argob('--collections', type=OpenFilePath, nargs='+',
                help="One or more, plottable data '.chp' files")
        self._create_argob('--plot_filename', type=SaveFilePath,
                help='Name of final plot file (valid extensions are '+\
                ".pdf, .png, .jpg)")
        self._create_argob('-e', '--expression_data', type=OpenFilePath,
                help='Path to the expression/gene data file, normally '+\
                "with '.exp' file extension. Must be tab delimited.")
        self._create_argob('--allow_probeset_many_gene', action='store_true',
                default=False, help='Allow probesets that map to multiple '+\
                'genes. Only used if probeset ids are present')
        self._create_argob('--gene_id_heading', default='gene',
                help='Column header identifier for the column containing '+\
                'Ensembl gene stable IDs')
        self._create_argob('--probeset_heading', default='probeset',
                help='Column header identifier for the column containing '+\
                'probeset IDs')
        self._create_argob('--expression_heading', default='exp',
                help='Column header identifier for the column containing '+\
                'expression or expression difference scores')
        self._create_argob('--significance_heading', default='sig',
                help='Column header identifier for the column containing '+\
                'significance values')
        self._create_argob('--p_value_heading', default='p_val',
                help='Column header identifier for the column containing '+\
                'expression difference p values')
        self._create_argob('--sep', default='\t', help='data field delimiter')

        self._create_argob('-s', '--sample', choices=\
                db_query.get_all_sample_names(self.db_path),
                help='Select an existing sample by name from the DB')
        self._create_argob('--samples', nargs='+', choices=\
                db_query.get_all_sample_names(self.db_path),
                help='Select one or more existing samples by '+\
                'space-separated name from the DB')
        # absolute or differential expression
        self._create_argob('--expr_sample', choices=\
                db_query.get_expr_diff_sample_names(self.db_path),
                help='Select an existing sample (absolute or differential '+\
                'expression) by name from the DB')
        # absolute expression only
        self._create_argob('--abs_expr_sample', choices=\
                db_query.get_expr_sample_names(self.db_path),
                help='Select an absolute expression sample by name from the DB')
        self._create_argob('--abs_expr_samples', choices=\
                db_query.get_expr_sample_names(self.db_path), nargs='+',
                help='Select one or more absolute expression samples by '+\
                'space-separated name from the DB')
        # differential expression samples only
        self._create_argob('--diff_sample', choices=\
                db_query.get_diff_sample_names(self.db_path),
                help='Select a differential expression sample by '+\
                'name from the DB')

        # args for creating new sample entries
        samplesStr = ', '.join(db_query.get_all_sample_names(self.db_path))
        self._create_argob('-n','--name', default=None, help='Enter a new '+\
                'sample name. Existing samples: ' + samplesStr)
        self._create_argob('-d', '--description', help='give your sample a '+\
                'description')
        self._create_argob('--sample_type',
                choices=sample_type_names, default='abs_expr',
                help='One of: '+ ', '.join(sample_type_desc)+\
                'Select the type of data you want entered from '+\
                ' '.join(sample_type_names))

        self._create_argob('--reffile1', type=OpenFilePath,
                help='Related file 1. First absolute expression sample in '+\
                'differential expression experiment. Full path to file.')
        self._create_argob('--reffile2', type=OpenFilePath,
                help='Related file 2. Second absolute expression sample in '+\
                'differential expression experiment. Full path to file.')

        # Export Centred Counts args
        self._create_argob('-B', '--BAMorBED', type=OpenFilePath,
                help='Path to counts data. Indexed BAM, BED, BEDgraph, '+\
                'or WIG file')

        self._create_argob('-f', '--overwrite', action='store_true',
                help='Overwrite any existing file(s)', default=False)
        self._create_argob('--tab_delimited', action='store_true',
                help='output to tab delimited format', default=False)
        self._create_argob('--chr_prefix', default='', help='String '+\
                "prefixing chromosome numbers/names. Often 'chr' or 'chr_'")

        self._create_argob('--wig', type=OpenFilePath,
                help='WIG file containing expression data to open')
        self._create_argob('--exp', type=SaveFilePath,
                help='.exp gene expression file save path')

        # Used to by plot_counts.py to save gene list
        self._create_argob('--write_genes_by_rank', type=SaveFilePath,
                help='Path to write the ENSEMBL gene ids seen in the '+\
                'plot, from highest expression rank to lowest. Will be '+\
                'modified by study names. A separate file per study will ' +\
                'be created.')

    def _add_sampling_args(self):
        """ All arguments relate to conditional selection of data """
        # chrom choice
        self._create_argob('-C', '--chrom', default=None,
                help='Restrict study to specific chromosome',
                choices=(db_query.get_chroms(self.db_path)) )

        # group genes into sets ranked by expression
        self._create_argob('-g', '--group_size', default='All',
                help='Number of genes to group together in one plotted '+\
                "line. Can be an integer number or 'All'", important=True)

        self._create_argob('--num_genes', type=int,
                help='Use this number of top-ranked genes from study')

        self._create_argob('--group_location', default='all',
                choices=['all', 'top', 'middle', 'bottom'],
                help='Show only a representative group of genes in a '+\
                'study based on their expression')

        self._create_argob('--ranks', action='store_true',
                help='Use rank-based data values instead of counts or '+\
                'expression')

        self._create_argob('--sample_extremes', type=float,
                default=0.0, help='The proportion of least- and most- '+\
                'extreme-value absolute expressed genes to treat '+\
                'separately. Set to 0.0 to disable.')

        # Use moving-mean (boxcar) smoothing on lines
        self._create_argob('--smoothing', type=int, default=0,
                help='Window size for smoothing of plot data in bases '+\
                '(integer), e.g. 50 is suitable for nucleosomes.')

        self._create_argob('--binning', type=int, default=0,
                help='Display counts-per-base within integer-sized bins '+\
                'of bases across the plot')

        # Filter out over- and under-expression outliers
        self._create_argob('--data_cutoff', type=float, default=0.0,
                help='Probability cutoff using one-sided Chebyshev filtering '+\
                'of extreme gene counts. Exclude genes if the probability '
                'of the observed tag count is less than or equal to this '+\
                'value e.g. 0.05. This implies values more than 4.36 standard '+\
                'deviations from the mean may be outliers with p < 0.5. p=0.01 '+\
                'gives a cut-off of 9.95 standard deviations.')
        # Filter gene counts in each line
        self._create_argob('--line_cutoff', type=float, default=0.0, help='Use '+\
                'Probability cutoff using two-sided Chebyshev filtering '+\
                'within each plottable group of genes, to remove outliers. '+\
                'Suggested cut-off p value = 0.01 implies genes with with '+\
                'values more than 7.07 standard deviations from the mean '+\
                'may be outliers.')

        # Export Centred Counts args
        self._create_argob('--feature_type', choices=['TSS', 'UTR_Exon',
                'Exon_Intron', 'Intron_Exon', 'Exon_UTR', 'Gene_3p'],
                help='Gene feature options: TSS, UTR_Exon, Exon_Intron, '+\
                'Intron_Exon, Exon_UTR, Gene_3p')

        self._create_argob('--window_upstream', type=int, default=1000,
                help='Region start distance relative to feature in '+\
                'bases, e.g. 1000')
        self._create_argob('--window_downstream', type=int, default=1000,
                help='Region finish distance relative to feature in '+\
                'bases, e.g. 1000')

        self._create_argob('--multitest_signif_val', type=int,
                help='Restrict plot to genes that pass multi-test '+\
                'significance. Valid values: 1, 0, -1', default=None)

        self._create_argob('--include_targets', default=None, nargs='+',
                help='Restrict the output to at most include the genes '+\
                'in the chosen Target Gene lists.', choices=\
                [str(s) for s in db_query.get_target_gene_names(self.db_path)])
        self._create_argob('--exclude_targets', default=None, nargs='+',
                help='Restrict the output to not includes any genes '+\
                'from the chosen Target Gene lists.', choices=\
                [str(s) for s in db_query.get_target_gene_names(self.db_path)])

    def _add_plot_args(self):
        """ Arguments specifically related to showing graphical plots """
        self._create_argob('-y', '--ylim', default=None,
                help='comma separated minimum-maximum yaxis values (e.g. 0,3.5)')
        self._create_argob('-H', '--fig_height', type=float, default=6.5*2.5,
                help='Figure height (cm)')
        self._create_argob('-W', '--fig_width', type=float, default=10*2.5,
                help='Figure width (cm)')

        # Important note, grid_lines are an absolute scale!
        self._create_argob('--xgrid_lines', type=float, default = 100,
                help='major grid-line spacing on x-axis')
        self._create_argob('--ygrid_lines', type=float, default = None,
                help='major grid-line spacing on y-axis')
        self._create_argob('--grid_off', action='store_true',
                help='Turn grid lines off')
        # tick spacing along axes
        self._create_argob('--xtick_interval', type=int, default=2,
                help='number of blank ticks between labels')
        self._create_argob('--ytick_interval', type=int, default=2,
                help='number of blank ticks between labels')
        # offsets ticks from edge of plot
        self._create_argob('--offset_ticks', action='store_true',
                help='offset ticks from edge of plot')
        # Smooth top and right borders for cleaner looking plot
        self._create_argob('--clean_plot', action='store_true', default=False,
                help='Remove tick marks and top and right borders ')
        # background colour - black or white
        self._create_argob('-b', '--bgcolor', default='black',
                help='Plot background color',
                choices=['black', 'white'])
        # side colour bar for 3rd-dimension range
        self._create_argob('--colorbar', action='store_true',
                help="Add colorbar to figure", default=False)

        # Create options for font sizes, colours and labels
        # Headline of the plot
        self._create_argob('--font', help='font name, defaults to Vera Sans')
        self._create_argob('--title', help='Plot title')
        self._create_argob('--title_size', type=int, default=18,
                help='font size for title')
        # Need options for size and style
        # Axis labels and font sizes
        self._create_argob('--ylabel', default = '',
             help='Label for the y-axis')
        self._create_argob('--xlabel', default = '',
                help='Label for the x-axis')
        self._create_argob('--xfont_size', type=int, default=14,
                help='font size for x label')
        self._create_argob('--yfont_size', type=int, default=14,
                help='font size for y label')
        # Optional axis units text
        self._create_argob('--yaxis_units', default='',
                help='Text showing units of y-axis of plot')
        self._create_argob('--xaxis_units', default='',
                help='Text showing units of x-axis of plot')
        self._create_argob('--xaxis_text1', default='',
                help='X-axis label for plot 1')
        self._create_argob('--xaxis_text2', default='',
                help='X-axis label for plot 2')

        # Turn on line legend and set font size
        self._create_argob('-l', '--legend', action='store_true', default=False,
                help='Automatically generate a figure legend.')
        self._create_argob('--legend_font_size', type=int, default=12,
                help='Point size for legend characters')

        # Vertical line showing the centred feature
        self._create_argob('--vline_style',
                default = '-.', choices=['-.', '-', '.'],
                help='line style for centred vertical line')
        self._create_argob('--vline_width', type=int,
                default = 2,
                help='line width for centred vertical line')

        # Plotted line opacity
        self._create_argob('--line_alpha', type=float, default=1.0,
                help='Opacity of lines')
        # Line width
        self._create_argob('--line_width', type=float, default=2.0,
                help='Thickness of plotted lines')

        # Plot as grey-scale
        self._create_argob('--grey_scale', action='store_true',
                help='Plot colour range as grey-scale')
        # Line colour range restriction
        self._create_argob('--restrict_colors', default='0.1,0.9',
                help='Comma separated 0.0..1.0 range for color '+\
                     'spectrum limits (e.g. 0.1,0.9)')

        # Plot as series of images
        self._create_argob('--plot_series', action='store_true',
                help='Plot series of figures. A directory called '+\
                'plot_filename-series will be created.')

        # Position of legend?
        self._create_argob('--text_coords', default=None,
                help='x, y coordinates of series text (e.g. 600,3.0)')

        self._create_argob('--div', help='Path to the plottable data, to divide '+\
                'other loaded plottable data')
        self._create_argob('--div_by', choices=['average', 'median', 'top'],
                help='divide by a single line rather than line-for-line')

        self._create_argob('--confidence_intervals', action='store_true',
                help='Show confidence intervals around plot lines')

        # Uses the 'tag count' arg value in the export file
        self._create_argob('--no_normalise', action='store_true',
                help="Don't normalise by Reads Per mapped-Million "+\
                "with 'mean counts' metric")

    def _add_db_args(self):
        """ options for starting a ChipPy DB """
        self._create_argob('--save_db_dir', type=DirPath,
                help='path to directory where chippy.db will be saved')
        self._create_argob('--save_db_prefix', default = '',
                help='Prefix string to DB file name')
        self._create_argob('--ensembl_release', type=int,
                help='Ensembl release to use.')
        self._create_argob('--species', default='mouse',
                help="Create for species e.g. 'mouse','human'")
        self._create_argob('--hostname', default=None,
                help='hostname for SQL Ensembl server')
        self._create_argob('--username', default=None,
                help='username SQL Ensembl server')
        self._create_argob('--password', default=None,
                help='password for SQL Ensembl server')
        self._create_argob('--port', default=None, type=int,
                help='Port for SQL Ensembl server')
        self._create_argob('--dummy_data', default=False, action='store_true',
                help='Create dummy expression data')

    def _add_misc_args(self):
        """ various options that don't fall into a category above """

        self._create_argob('-m', '--counts_metric',
                choices=['mean', 'frequency', 'stdev'], default='mean',
                help='Select the metric for combining counts. Mean '+\
                'is RPM normalised by default. Frequency is Counts per '+\
                'Million Base Counts window-normalised and may be '+\
                'helpful  when reads are expected to have mis-mapped '+\
                'due to biological confounding. Stdev plots in standard '+\
                'deviations from the mean count and is useful for showing '+\
                'trends in feature exaggeration or positional movement.')

        self._create_argob('-t', '--test_run', action='store_true',
                help="Test run, don't write output",
                default=False)

        # should never be displayed with the app
        self._create_argob('--show_log', action='store_true',
                help='Print recent log entries following execution',
                display=False)

        self._create_argob('--max_chrom_size', type=int, default=300000000,
                help='Largest possible chromosome size in bases')

        # Used by some auxiliary scripts
        self._create_argob('--x_axis_is_log', action='store_true',
            help='Plot x-axis with log2 values')
        self._create_argob('--y_axis_is_log', action='store_true',
            help='Plot y-axis with log2 values')

    def _add_diff_abs_plots_specific_args(self):
        """ These args are specifically for the diff_abs_plots script """
        # absolute expression only
        self._create_argob('--abs_expr_sample1', choices=\
        db_query.get_expr_sample_names(self.db_path),
            help='Select an absolute expression sample to use')
        self._create_argob('--abs_expr_sample2', choices=\
        db_query.get_expr_sample_names(self.db_path),
            help='Select an absolute expression sample to use')

        # dot colouring
        self._create_argob('--extremes_colour', default='blue', choices=['blue',
                'red', 'yellow', 'green', 'magenta', 'orange', 'cyan'],
                help='Colour of dots for absolute expression marked '\
                'as extreme.')
        self._create_argob('--signif_colour', default='blue', choices=['blue',
                'red', 'yellow', 'green', 'magenta', 'orange', 'cyan'],
                help='Colour of dots for difference of expression marked '\
                'as significant.')
        self._create_argob('--bulk_colour', default='blue', choices=['blue',
                'red', 'yellow', 'green', 'magenta', 'orange', 'cyan'],
                help='Colour of dots for all relatively unexceptional '\
                'expression values.')

        # hide unwanted plot areas
        self._create_argob('--hide_extremes', default=False, action='store_true',
                help='Do not show absolute expression considered extreme')
        self._create_argob('--hide_signif', default=False, action='store_true',
                help='Do not show difference expression considered '+\
                'significant')
        self._create_argob('--hide_bulk', default=False, action='store_true',
                help='Do not show expression values considered normal')

        # misc. options that don't fit other categories
        self._create_argob('--plot1_name',
                default=None, help='Output path for first plot')
        self._create_argob('--plot2_name',
                default=None, help='Output path for second plot')

    def _add_counts_vs_expr_specific_args(self):
        """ These args are specifically for the counts_vs_expr script """
        self._create_argob('--x_axis_type', choices=['expression', 'counts'],
                help='Select the data unit type for the x-axis',
                default='expression')

        self._create_argob('--region_feature', choices=['total', 'promoter',
                'coding', 'feature'], help='Which part of a counts region '+\
                'should be used to calculate rank and score',
                default='total')

        self._create_argob('--counts_is_ranks', action='store_true',
                help='Plot chromatin counts as ranks rather than absolute values')
        self._create_argob('--expr_is_ranks', action='store_true',
                help='Plot expression as ranks rather than absolute values')

    def _add_counts_distribution_specific_args(self):
        """ Args specific to the counts_distribution script """
        # Also used by expr_distribution
        self._create_argob('--plot_type', choices=['dot', 'line', 'box',
                'hist'], default='dot',
                help='Type of plot for comparing counts data')
        self._create_argob('--counts_region', choices=['feature', 'total',
                'promoter', 'coding'], default='feature', help=\
                'Region of export window for calculating score')

    def _add_expr_distribution_specific_args(self):
        """ Args specific to the expr_distribution script """
        pass
        # would also use "plot_type"

    def getReqCogentOpts(self):
        return [argob.asCogentOpt() for argob in self.argobs\
                if argob.required]

    def getOptCogentOpts(self):
        return [argob.asCogentOpt() for argob in self.argobs\
                if not argob.required]

    def __init__(self, positional_args=None, required_args=None,
            optional_args=None):
        """ calls _create_argobs on every possible argument """
        self.parser = argparse.ArgumentParser(version='ChipPy r'+str(__version__))
        self.required_args = required_args
        self.optional_args = optional_args

        self.argobs = []
        self.db_path = None

        # We need to handle the 'db_path' positional argument ourselves
        # as several arguments require it already loaded
        if positional_args:
            if 'db_path' in positional_args:
                # create a positional, non-displayed ArgOb
                for arg in sys.argv:
                    if not arg.startswith('-'):
                        # it's positional so take this to be db_path
                        possible_db_path = str(arg).strip()
                        if 'chippy' in possible_db_path.lower() and\
                                '.db' in possible_db_path.lower():
                            self.db_path = possible_db_path
                            db = ArgOb('db_path', type=OpenFilePath,
                                    help='Path to ChippyDB', display=False,
                                    default=possible_db_path)
                            self.argobs.append(db)
                            self.required_args.append(db)
                            print self.db_path, 'selected as ChipPy database'
                            break

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

        # args specific to the individual scripts
        self._add_diff_abs_plots_specific_args()
        self._add_counts_vs_expr_specific_args()
        self._add_counts_distribution_specific_args()
        self._add_expr_distribution_specific_args()
