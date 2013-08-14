from __future__ import division
import warnings
import sys
sys.path.extend(['..', '../..'])
import numpy
from matplotlib import pyplot, rc, cm, font_manager
from matplotlib.mpl import colorbar
from matplotlib.ticker import MultipleLocator

from cogent.util.progress_display import display_wrap
from chippy.util.run_record import RunRecord
from math import log10, floor, ceil
from numpy import NINF, PINF

ColorbarBase = colorbar.ColorbarBase

__author__ = 'Gavin Huttley, Cameron Jack'
__copyright__ = 'Copyright 2011-2013, Gavin Huttley, Cameron Jack, Anuj Pahwa'
__credits__ = ['Gavin Huttley', 'Cameron Jack']
__license__ = 'GPL'
__maintainer__ = 'Cameron Jack'
__email__ = 'cameron.jack@anu.edu.au'
__status__ = 'pre-release'
__version__ = '0.1'

class _Plottable(object):
    """ Base class for handling plotting. Defines the appearance
        of a plot.
    """

    def __init__(self, height, width, bgcolor, grid_off, pad=10,
            xaxis_lims=None, yaxis_lims=None, xy_tick_spaces=None,
            xy_tick_intervals=None, linewidth=2,
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

        self.fig = None
        self.ax = None
        self._legend_patches = []
        self._legend_labels = []
        self._line_collection = []
        self._colorbar = colorbar
        self.clean = clean

    ### private helper methods

    def _auto_grid_lines(self, y_ceiling, test_run=False):
        """ Returns a float that is a 'round' looking number to use for
            the grid lines
        """
        rr = RunRecord('_auto_grid_lines')
        if y_ceiling > 0:
            ypower = log10(y_ceiling)
            if ypower < 0:
                rounding_places = 0 - int(floor(ypower))
                y_ceiling = float(ceil(y_ceiling*(10**rounding_places))/\
                                  (10**rounding_places))
                grid_line_val = y_ceiling/10.0
            else:
                y_ceiling = ceil(y_ceiling)
                if y_ceiling <= 10:
                    grid_line_val = round(y_ceiling/10.0, 1)
                else:
                    grid_line_val = y_ceiling/10.0
        else:
            rr.dieOnCritical('Inappropriate y-axis value', y_ceiling)
        if test_run:
            rr.addInfo('Y-grid-line spacing', '%e' % grid_line_val)
        return grid_line_val

    def _auto_y_lims(self, y=None, plot_lines=None, rounding=True,
            test_run=False):
        """ Takes either a list of y values or a list of plotlines.
            Returns ylims(y_min_limit, y_max_limit)
            Cannot have negative y-axis.
            Defaults to min = 0.0, max = 1.0
        """
        rr = RunRecord('_auto_y_lims')
        # Get min/max y-axis values
        y_min_limit = PINF
        y_max_limit = NINF
        if plot_lines:
            for line in plot_lines:
                y_min_limit = min(min(line.counts), y_min_limit)
                y_max_limit = max(max(line.counts), y_max_limit)
        elif y:
            if len(y) > 1:
                for counts_array in y:
                    y_min_limit = min(min(counts_array), y_min_limit)
                    y_max_limit = max(max(counts_array), y_max_limit)
            else: # just a single counts array, so single line
                y_min_limit = min(y[0])
                y_max_limit = max(y[0])
        else:
            rr.dieOnCritical('No y-array or plotlines provided', 'Failed')

        if rounding:
            # Round min/max values to whole values for nice plots

            # For fractional counts then scale the rounding appropriately
            if y_max_limit > 0:
                ypower = log10(y_max_limit) # check scale

                if ypower < 0:
                    rounding_places = 0 - int(floor(ypower))
                    y_ceiling = float(ceil(y_max_limit * (10**rounding_places))/\
                                  (10**rounding_places))
                    y_floor = float(floor(y_min_limit * (10**rounding_places))/\
                                (10**rounding_places))
                elif ypower == 0:
                    y_floor = 0.0
                    y_ceiling = 1.0
                else:
                    # round up to 2 significant digits
                    ypower = ceil(log10(y_max_limit))
                    y_ceiling = ceil( y_max_limit/(10**(ypower-1)) ) * (10**(ypower-1))
                    y_floor = floor(y_min_limit)
            elif y_max_limit == 0:
                y_floor = 0.0
                y_ceiling = 1.0
            else:
                rr.dieOnCritical('Negative max y-axis value', y_max_limit)
        else:
            y_floor = y_min_limit
            y_ceiling = y_max_limit

        if test_run:
            rr.addInfo('Y-axis min', y_min_limit)
            rr.addInfo('Y-axis max', y_max_limit)
            rr.addInfo('Y-axis auto floor', y_floor)
            rr.addInfo('Y-axis auto ceiling', y_ceiling)

        return (y_floor, y_ceiling)

    ### private methods called internally as needed
    
    def _get_figure_and_axes(self, title=None, xlabel=None, ylabel=None):
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
                d, r = divmod(i, self.xtick_interval)
                if r != 0:
                    xtick.set_visible(False)
        
        if self.ytick_interval is not None:
            yticks = ax.yaxis.get_major_ticks()
            for i, ytick in enumerate(yticks):
                d, r = divmod(i, self.ytick_interval)
                if r != 0:
                    ytick.set_visible(False)
        
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

        Sets the background to either black or white. White plots are
        designed for a minimal, 'clean' look.
        bgcolor = 'black'|'white'
        vline = (width, style, color)
        """
        x, vline_width, vline_style, vline_color = vline
        if bgcolor.lower() == 'black':
            if grid_off is True:
                self.grid = False
                vline_color = 'k'
            else:
                self.grid = {'color': 'w'}
                vline_color = 'w'
            self.bgcolor='0.1'
        else:
            if grid_off is True:
                self.grid = False
                vline_color = 'w'
            else:
                self.grid = {'color': 'k'}
                vline_color = 'k'
            self.bgcolor = '1.0'

        self.vline = dict(x=x, linewidth=vline_width,
                linestyle=vline_style, color=vline_color)

    def _set_axes(self, y_vals=None, plot_lines=None, test_run=False):
        """ Gets called by the __call__ method but is also available for
        re-scaling of plots.

        1) Set the axes to y_min_limit and y_max_limit or call
        auto-calculate.
        2) Set y-tick-space or auto-calculate
        """
        rr = RunRecord('_set_axes')
        if self.ylims is None:
            if plot_lines or y_vals: # auto-calculate y-axis min and max
                self.ylims = self._auto_y_lims(y=y_vals,
                        plot_lines=plot_lines, test_run=test_run)
            else:
                rr.dieOnCritical('y data or plotlines', 'Values missing')

        y_min_limit, y_max_limit = self.ylims
        # set grid-lines/tick marks
        if self.ytick_space is None:
            self.ytick_space = self._auto_grid_lines(y_max_limit, test_run=test_run)

        if self.ytick_interval is None:
            # If self.ytick_space is even, then set to 2, otherwise 1.
            if self.ytick_space%2 == 0:
                self.ytick_interval = 2
            else:
                self.ytick_interval = 1

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

    def check_y_axis_scale(self, maxY=None, plot_lines=None):
        rr = RunRecord('check_y_axis_scale')
        if plot_lines:
            maxY = 0
            for line in plot_lines:
                maxY = max(line.getMaxCount(), maxY)

        if self.ylims is not None:
            if maxY > self.ylims[1]:
                rr.addWarning('ylimit may be too small, ymax=', str(maxY))
            elif maxY*2 < self.ylims[1]:
                rr.addWarning('ylimit may be too large, ymax=', str(maxY))
        else:
            rr.addWarning('y-axis limits', 'Not set')

### Public classes implementing Plottable

class PlottableSingle(_Plottable):
    """Plots a single line"""
    def __init__(self, *args, **kwargs):
        super(PlottableSingle, self).__init__(*args, **kwargs)
    
    def __call__(self, x, y, plot_lines=None, stderr=None, color=None,
            cmap='RdBu', clean=False, xlabel=None, ylabel=None, title=None,
            label=None):
        rr = RunRecord('PlottableSingle__call__')
        cmap = getattr(cm, cmap)
        if color is None:
            color = 'b'
        elif type(color) != str:
            color = cmap(color)

        self._set_axes(y_vals=y, plotlines=None, test_run=False)

        self.fig, self.ax = self._get_figure_and_axes(title=title,
                xlabel=xlabel, ylabel=ylabel)

        self.clean=clean

        self.check_y_axis_scale(y)
            
        patches = pyplot.plot(x, y, color=color, linewidth=self.linewidth,
                label=label)
        self._legend_patches.append(patches[0])
        self._legend_labels.append(label)
        
        if stderr is not None:
            upper = 1.96 * stderr + y
            lower = -1.96 * stderr + y
            pyplot.fill_between(x, upper, lower, alpha=0.2, color=color)

class PlottableGroups(_Plottable):
    """plot groups of data on the same panel"""
    def __init__(self, *args, **kwargs):
        super(PlottableGroups, self).__init__(*args, **kwargs)
    
    @display_wrap
    def __call__(self, x, y_series=None, plot_lines=None, color_series=None,
            alpha=None, series_labels=None, label_coords=None, cmap=None,
            colorbar=False, clean=False, xlabel=None, ylabel=None,
            title=None, filename_series=None, labels=None, labels_size=None,
            stderr=None, show_legend=False, plot_CI=False, ui=None):
        rr = RunRecord('PlottableGroups__call__')

        if not y_series and not plot_lines:
            rr.dieOnCritical('No data supplied', 'Failed')

        bbox = dict(facecolor='b', alpha=0.5)

        if plot_lines:
            # Get all data from plot_lines
            plot_lines = sorted(plot_lines, key=lambda line: line.rank, reverse=True)
            y_series = [line.counts for line in plot_lines]
            color_series = [line.color for line in plot_lines]
            labels = [line.label for line in plot_lines]

            # Reverse ranks so that rank 0 is colored red, must be in range 0..1
            ranks = sorted([line.rank/len(plot_lines) for line in plot_lines], reverse=True)
        if cmap:
            cmap_r = getattr(cm, '%s_r' % cmap)
            cmap = getattr(cm, cmap)

        # All data comes from separate arrays
        if color_series is not None:
            assert len(y_series) == len(color_series)
        
        if series_labels is not None:
            assert len(y_series) == len(series_labels)
            assert label_coords is not None
            label_x, label_y = label_coords

        self._set_axes(y_vals=y_series, plot_lines=plot_lines, test_run=False)

        self.fig, self.ax = self._get_figure_and_axes(title=title,
                xlabel=xlabel, ylabel=ylabel)

        self.clean=clean

        if self._colorbar and colorbar:
            # probably need to set a limit on how big this will be
            ax2 = self.fig.add_axes([0.925, 0.1, 0.025, 0.8])
            cb = ColorbarBase(ax2, ticks=[0.0, 1.0], cmap=cmap_r,
                    orientation='vertical')
            cb.set_ticklabels(['Low', 'High'])
            ax = self.fig.sca(self.ax) # need to make main axis the current axis again

        num = len(y_series)

        for i in ui.series(range(num), noun='Applying lines to plot'):
            if color_series is None:
                color = 'b'
            elif cmap and plot_lines:
                color = cmap(ranks[i])
            elif type(color_series[i]) != str:
                color = cmap(color_series[i])
            elif type(color_series[i]) == str:
                color = color_series[i]
            
            if series_labels is not None:
                txt = ax.text(label_x, label_y, series_labels[i],
                        bbox=bbox, color='w', fontsize=self.legend_font_size)

            if color_series is not None:
                y = y_series[i]
            else:
                y = y_series

            if labels is not None and show_legend:
                self._legend_labels.append(labels[i])
                patches, = pyplot.plot(x, y, color=color,
                        linewidth=self.linewidth, label=labels[i])
                self._legend_patches.append(patches)
                self.legend(labels_size)
            else:
                pyplot.plot(x, y, color=color, linewidth=self.linewidth, alpha=alpha)

            if filename_series is not None:
                pyplot.savefig(filename_series[i])
            
            if series_labels is not None:
                ax.texts.remove(txt)

            # Show confidence interval around each line
            if plot_CI and plot_lines is not None:
                upper = 1.96 * plot_lines[i].stderr + y_series[i]
                lower = -1.96 * plot_lines[i].stderr + y_series[i]
                pyplot.fill_between(x, upper, lower, alpha=0.3, color='green')

        self.check_y_axis_scale(maxY=max(y), plot_lines=plot_lines)

