from __future__ import division
import warnings
import sys
sys.path.extend(['..', '../..'])
with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    from matplotlib import pyplot, rc, cm, font_manager
    from matplotlib.mpl import colorbar
    from matplotlib.ticker import MultipleLocator
    from cogent.util.progress_display import display_wrap
from chippy.util.run_record import RunRecord
from math import log10, floor, ceil

ColorbarBase = colorbar.ColorbarBase

__author__ = 'Gavin Huttley, Cameron Jack'
__copyright__ = 'Copyright 2011-2013, Gavin Huttley, Cameron Jack, Anuj Pahwa'
__credits__ = ['Gavin Huttley', 'Cameron Jack']
__license__ = 'GPL'
__maintainer__ = 'Cameron Jack'
__email__ = 'cameron.jack@anu.edu.au'
__status__ = 'pre-release'
__version__ = '0.2'

class FigureDetails(object):
    """
        A 'lite' Plottable object to aid in passing useful information to
        plotting code. Should likely be merged with _Plottable.
    """

    def __init__(self, x_size=5, y_size=3, title=None, x_text=None,
                 y_text=None):
        self.x_size = x_size
        self.y_size = y_size
        self.title = title
        self.x_text = x_text
        self.y_text = y_text


class _Plottable(object):
    """ Base class for handling plotting. Defines the appearance
        of a plot.
    """

    def __init__(self, height, width, bgcolor, grid_off, pad=10,
            xaxis_lims=None, yaxis_lims=None, xy_tick_spaces=None,
            xy_tick_intervals=None, offset_ticks=False, linewidth=2,
            xy_label_fontsizes=(12,12), vline=None,
            legend_font_size=10,
            ioff=None, colorbar=False, clean=False):
        """
        height, width = physical size of plot in inches
        bgcolor = background color {black | white}
        grid_off = True|False (default False)
        pad = tick mark padding
        xaxis_lims = (x_min, x_max)
        yaxis_lims = (y_min, y_max)
        xy_tick_spaces = (x, y) tick spacing
        xy_tick_intervals = (x, y) display values for ticks every n (int)
        linewidth = thickness of plot lines
        xy_label_fontsizes = (x, y) font size for axis labels
        vline = (x, width, style, color)
        legend_font_size = font size for the plot legend
        ioff = interactive plot (True is passed in by default)
        colorbar = include a color scale bar with the plot
        clean = removes top and right plot edges and their tick marks
        """
        super(_Plottable, self).__init__()
        if ioff is not None:
            pyplot.ioff()

        rc('xtick.major', pad=pad)
        rc('xtick.minor', pad=pad)
        
        self.height = height
        self.width = width
        self._set_background(bgcolor, grid_off, vline)

        self.xlims = xaxis_lims
        self.ylims = yaxis_lims

        self.vline = vline
        self.legend_font_size = legend_font_size
        self.linewidth = linewidth

        self.xlabel_fontsize, self.ylabel_fontsize = xy_label_fontsizes
        self.xtick_space, self.ytick_space = xy_tick_spaces
        self.xtick_interval, self.ytick_interval = xy_tick_intervals
        self.offset_ticks = offset_ticks

        self.fig = None
        self.ax = None
        self._legend_patches = []
        self._legend_labels = []
        self._line_collection = []
        self._colorbar = colorbar
        self.clean = clean

    ### private helper methods

    def _auto_grid_lines(self, y_diff, test_run=False):
        """ Returns a float that is a 'round' looking number to use for
            the grid lines
        """
        rr = RunRecord('_auto_grid_lines')
        if y_diff > 0:
            ypower = log10(y_diff)
            if ypower < 0:
                rounding_places = 0 - int(floor(ypower))
                y_diff = float(ceil(y_diff*(10**rounding_places))/\
                                  (10**rounding_places))
                grid_line_val = y_diff/10.0
            else:
                y_ceiling = ceil(y_diff)
                if y_ceiling <= 10:
                    grid_line_val = round(y_ceiling/10.0, 1)
                else:
                    grid_line_val = y_ceiling/10.0
        else:
            rr.dieOnCritical('Y-axis length must be greater than 0', y_diff)
        if test_run:
            rr.addInfo('Y-grid-line spacing', '%e' % grid_line_val)
        return grid_line_val

    def _auto_y_lims(self, minY, maxY, rounding=True, test_run=False):
        """
            Takes a list of plotlines.
            Returns ylims(y_min_limit, y_max_limit)
            Defaults to min = 0.0, max = 1.0
        """
        rr = RunRecord('_auto_y_lims')

        y_floor = minY
        y_ceiling = maxY
        if rounding:
            # Round min/max values to whole values for nice plots

            # For fractional counts then scale the rounding appropriately
            if maxY > 0:
                ypower = log10(maxY) # check scale

                if ypower < 0:
                    rounding_places = 0 - int(floor(ypower))
                    y_ceiling = float(ceil(maxY * (10**rounding_places))/
                                  (10**rounding_places))
                    y_floor = float(floor(minY * (10**rounding_places))/
                                (10**rounding_places))
                elif ypower == 0:
                    y_floor = 0.0
                    y_ceiling = 1.0
                else:
                    # round up to 2 significant digits
                    ypower = ceil(log10(maxY))
                    y_ceiling = ceil( maxY/(10**(ypower-1)) ) * (10**(ypower-1))
                    y_floor = floor(minY)
            elif maxY == 0:
                y_floor = 0.0
                y_ceiling = 1.0
            else:
                rr.dieOnCritical('Negative max y-axis value', maxY)

        if test_run:
            rr.addInfo('Y-axis min', minY)
            rr.addInfo('Y-axis max', maxY)
            rr.addInfo('Y-axis auto floor', y_floor)
            rr.addInfo('Y-axis auto ceiling', y_ceiling)

        return tuple([y_floor, y_ceiling])

    def getFigureAndAxes(self, title=None, xlabel=None, ylabel=None):
        """returns the figure and axis ready for display"""
        if self.fig is not None:
            return self.fig, self.ax
        
        if self.xlabel_fontsize:
            rc('xtick', labelsize=self.xlabel_fontsize)
        if self.ylabel_fontsize:
            rc('ytick', labelsize=self.ylabel_fontsize)
        
        fig = pyplot.figure(figsize=(self.width, self.height))
        
        if self._colorbar:
            ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        else:
            ax = pyplot.gca()
        
        ax_kwargs = {}
        if self.xlims is not None:
            ax_kwargs['xlim'] = self.xlims
        
        if self.ylims is not None:
            ax_kwargs['ylim'] = self.ylims
            
        pyplot.setp(ax, **ax_kwargs)
        
        if self.xtick_space is not None:
            major_locator = MultipleLocator(self.xtick_space)
            ax.xaxis.set_major_locator(major_locator)
        
        if self.ytick_space is not None:
            major_locator = MultipleLocator(self.ytick_space)
            ax.yaxis.set_major_locator(major_locator)
        
        if self.bgcolor is not None:
            ax.set_axis_bgcolor(self.bgcolor)

        if self.xtick_interval is not None:
            xticks = ax.xaxis.get_major_ticks()
            for i, xtick in enumerate(xticks):
                if self.offset_ticks:
                    d, r = divmod(i-1+self.xtick_interval/2, self.xtick_interval)
                else:
                    d, r = divmod(i-1, self.xtick_interval)
                xtick.set_visible(False)
                if r == 0:
                    xtick.set_visible(True)
        
        if self.ytick_interval is not None:
            yticks = ax.yaxis.get_major_ticks()
            for i, ytick in enumerate(yticks):
                if self.offset_ticks:
                    d, r = divmod(i-1+self.ytick_interval/2, self.ytick_interval)
                else:
                    d, r = divmod(i-1, self.ytick_interval)
                ytick.set_visible(False)
                if r == 0:
                    ytick.set_visible(True)

        if self.vline is not None:
            # e.g. x=0, ymin=0, ymax=1, linewidth=3, linestyle='-.', color='w'
            ax.axvline(**self.vline)
        
        if self.grid:
            ax.grid(**self.grid)
        
        if title:
            pyplot.title(title)
        
        if ylabel:
            pyplot.ylabel(ylabel, fontsize=self.ylabel_fontsize+2)
        
        if xlabel:
            pyplot.xlabel(xlabel, fontsize=self.xlabel_fontsize+2)

        ax.ticklabel_format(scilimits=(-2,4), axis='y')
        ax.ticklabel_format(scilimits=(-5,5), axis='x')
        self.fig = fig

        if self.clean is True:
            for loc, spine in ax.spines.iteritems():
                if loc in ['right','top']:
                    spine.set_color('none') # don't draw spine
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')

        self.ax = ax

        return self.fig, self.ax

    def _set_background(self, bgcolor, grid_off, vline):
        """ Called during initialisation.

        Sets the background to either black or white.
        bgcolor = 'black'|'white'
        vline = (x, width, style, color)
        """
        x, vline_width, vline_style, vline_color = vline
        if bgcolor.lower() == 'black':
            if grid_off is True:
                self.grid = False
            else:
                self.grid = {'color': 'w'}
                vline_color = 'w'
            self.bgcolor='0.0'
        else:
            if grid_off is True:
                self.grid = False
            else:
                self.grid = {'color': 'k'}
                vline_color = 'k'
            self.bgcolor = '1.0'

        if not grid_off:
            self.vline = dict(x=x, linewidth=vline_width,
                    linestyle=vline_style, color=vline_color)

    def setAxes(self, plot_lines, plot_CI=False, test_run=False):
        """
            Gets called by the __call__ method but is also available for
            re-scaling of plots.

            1) Set the axes to y_min_limit and y_max_limit or call
                auto-calculate.
            2) Set y-tick-space or auto-calculate
        """
        rr = RunRecord('setAxes')

        if not self.ylims:
            minY = self.getMinY(plot_lines, plot_CI)
            maxY = self.getMaxY(plot_lines, plot_CI)
            self.ylims = self._auto_y_lims(minY, maxY, test_run=test_run)

        y_min_limit, y_max_limit = self.ylims
        # set grid-lines/tick marks
        if not self.ytick_space:
            self.ytick_space = self._auto_grid_lines(y_max_limit-y_min_limit,
                    test_run=test_run)

        if not self.ytick_interval:
            self.ytick_interval = 2

        rr.addInfo('Y-max plot limit', '{:e}'.format(y_max_limit))
        rr.addInfo('Y-min plot limit', '{:e}'.format(y_min_limit))
        rr.addInfo('Y-grid-line spacing', '{:e}'.format(self.ytick_space))

    ### public methods for detailing a plottable object

    def ion(self):
        pyplot.ion()
    
    def show(self):
        pyplot.show()
    
    def savefig(self, filename, image_format='pdf'):
        pyplot.savefig(filename, format=image_format)
    
    def legend(self, fontsize=None):
        if self._legend_patches:
            if fontsize is None:
                prop = font_manager.FontProperties(size=self.xlabel_fontsize)
            else:
                prop = font_manager.FontProperties(size=fontsize)
            pyplot.legend(self._legend_patches, self._legend_labels,
                    prop=prop)

    def checkYAxisScale(self, plot_lines, plot_CI=False):
        """ Compare the set y-axis limits to the actual limits of the data """
        rr = RunRecord('checkYAxisScale')

        maxY = self.getMaxY(plot_lines, plot_CI)
        minY = self.getMinY(plot_lines, plot_CI)

        if self.ylims is not None:
            if maxY > self.ylims[1]:
                rr.addWarning('ylimit may be too small, ymax=', str(maxY))
            elif maxY*2 < self.ylims[1]:
                rr.addWarning('ylimit may be too large, ymax=', str(maxY))

            if minY < self.ylims[0]:
                rr.addWarning('ylimit may be too small, ymin=', str(minY))
            elif minY/2 > self.ylims[0]:
                rr.addWarning('ylimit may be too large, ymin=', str(minY))
        else:
            rr.addWarning('y-axis limits', 'Not set')


    def getMaxY(self, plot_lines, plot_CI=False):
        maxY = 0 # plots are never totally negative
        for line in plot_lines:
            peak = line.getMaxCount(include_stderr=plot_CI, se_adjust=1.96)
            if peak > maxY:
                maxY = peak
        return maxY

    def getMinY(self, plot_lines, plot_CI=False):
        minY = plot_lines[0].counts[0] # better than starting at zero
        for line in plot_lines:
            peak = line.getMinCount(include_stderr=plot_CI, se_adjust=1.96)
            if peak < minY:
                minY = peak
        return minY

### Public classes implementing Plottable

class PlottableSingle(_Plottable):
    """Plots a single line"""
    def __init__(self, *args, **kwargs):
        super(PlottableSingle, self).__init__(*args, **kwargs)

    @display_wrap
    def __call__(self, x_array, plot_lines=None, clean=False, xlabel=None,
            ylabel=None, title=None, plot_CI=False, ui=None):
        rr = RunRecord('PlottableSingle__call__')

        self.setAxes(plot_lines, plot_CI=plot_CI, test_run=False)
        self.checkYAxisScale(plot_lines, plot_CI=plot_CI)

        self.fig, self.ax = self.getFigureAndAxes(title=title,
              xlabel=xlabel, ylabel=ylabel)

        self.clean=clean

        for i, line in ui.series(enumerate(sorted(plot_lines,
                key=lambda line: (line.study,line.rank), reverse=True)),
                noun='Applying lines to plot'):
            self.ax.plot(x_array, line.counts, color=line.color,
                    linewidth=self.linewidth)

            # Show confidence interval around each line
            if plot_CI:
                #set shading alpha
                alpha = line.color[3]
                if alpha is None:
                    alpha = 0.9
                upper = 1.96 * line.stderr + line.counts
                lower = -1.96 * line.stderr + line.counts
                self.ax.fill_between(x_array, upper, lower, alpha=alpha/2.5,
                        color=line.color)

class PlottableGroups(_Plottable):
    """plot groups of data on the same panel"""
    def __init__(self, *args, **kwargs):
        super(PlottableGroups, self).__init__(*args, **kwargs)

    @display_wrap
    def __call__(self, x_array, plot_lines,
            colorbar=False, clean=False, xlabel=None, ylabel=None,
            title=None, filename_series=None, labels_size=None,
            show_legend=False, plot_CI=False, ui=None):
        rr = RunRecord('PlottableGroups__call__')

        if not plot_lines:
            rr.dieOnCritical('No data supplied', 'Failed')

        self.setAxes(plot_lines, plot_CI=plot_CI, test_run=False)
        self.checkYAxisScale(plot_lines, plot_CI=plot_CI)

        self.fig, self.ax = self.getFigureAndAxes(title=title,
             xlabel=xlabel, ylabel=ylabel)

        self.clean=clean

        if colorbar:
            # probably need to set a limit on how big this will be
            ax2 = self.fig.add_axes([0.925, 0.1, 0.025, 0.8])
            cb = ColorbarBase(ax2, ticks=[0.0, 1.0], cmap=cm.RdBu,
                    orientation='vertical')
            cb.set_ticklabels(['Low', 'High'])
            ax = self.fig.sca(self.ax) # need to make main axis the current axis again

        legend_lines = {}
        for i, line in ui.series(enumerate(sorted(plot_lines,
                key=lambda line: (line.study,line.rank), reverse=True)),
                noun='Applying lines to plot'):
            self.ax.plot(x_array, line.counts, color=line.color,
                    linewidth=self.linewidth)

            if show_legend:
                if line.study in legend_lines.keys():
                    if line.rank < legend_lines[line.study].rank:
                        legend_lines[line.study] = line
                else:
                    legend_lines[line.study] = line

            #if filename_series is not None:
            #    pyplot.savefig(filename_series[i])

            # Show confidence interval around each line
            if plot_CI:
                #set shading alpha
                alpha = line.color[3]
                if alpha is None:
                    alpha = 0.9
                upper = 1.96 * line.stderr + line.counts
                lower = -1.96 * line.stderr + line.counts
                self.ax.fill_between(x_array, upper, lower, alpha=alpha/2.5,
                        color=line.color)

        if show_legend:
            self.legend(labels_size)

            l_lines = [line for line in sorted(legend_lines.values(),
                    key=lambda x: x.study)]
            p_lines = [self.ax.plot(x_array, l.counts, color=l.color,
                    linewidth=self.linewidth)[0] for l in l_lines]
            study_names = [l.study for l in l_lines]

            self.ax.legend(p_lines, study_names)