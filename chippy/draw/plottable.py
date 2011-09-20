from __future__ import division
import warnings

from matplotlib import pyplot, rc, cm, font_manager
from matplotlib.mpl import colorbar
from matplotlib.ticker import MultipleLocator


from cogent.util.progress_display import display_wrap

ColorbarBase = colorbar.ColorbarBase

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"
__version__ = '0.1'

class _Plottable(object):
    """base class for handling plotting"""
    def __init__(self, height, width, bgcolor, grid, pad=10, ylim=None,
                xlim=None, xtick_space=None, ytick_space=None,
                xtick_interval=None, ytick_interval=None, linewidth=2,
                xlabel_fontsize=None, ylabel_fontsize=None, vline=None,
                ioff=None, colorbar=False):
        super(_Plottable, self).__init__()
        if ioff is not None:
            pyplot.ioff()
        
        rc('xtick.major', pad=pad)
        rc('xtick.minor', pad=pad)
        
        self.height = height
        self.width = width
        self.bgcolor = bgcolor
        self.grid = grid
        
        self.xlim = xlim
        self.ylim = ylim
        self.vline = vline
        
        self.xlabel_fontsize = xlabel_fontsize
        self.ylabel_fontsize = ylabel_fontsize
        
        self.xtick_space = xtick_space
        self.ytick_space = ytick_space
        self.xtick_interval = xtick_interval
        self.ytick_interval = ytick_interval
        self.linewidth = linewidth
        self.fig = None
        self.ax = None
        self._legend_patches = []
        self._legend_labels = []
        self._line_collection = []
        self._colorbar = colorbar
    
    def _get_figure_axis(self, title=None, xlabel=None, ylabel=None):
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
        if self.xlim is not None:
            ax_kwargs['xlim'] = self.xlim
        
        if self.ylim is not None:
            ax_kwargs['ylim'] = self.ylim
            
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
        
        ax.ticklabel_format(scilimits=(-3,4))
        self.fig = fig
        self.ax = ax
        return self.fig, self.ax
    
    def ion(self):
        pyplot.ion()
    
    def show(self):
        pyplot.show()
    
    def savefig(self, filename):
        pyplot.savefig(filename)
    
    def legend(self):
        if self._legend_patches:
            prop = font_manager.FontProperties(size=self.xlabel_fontsize)
            pyplot.legend(self._legend_patches, self._legend_labels,
                            prop=prop)
    

class PlottableSingle(_Plottable):
    """Plots a single line"""
    def __init__(self, *args, **kwargs):
        super(PlottableSingle, self).__init__(*args, **kwargs)
    
    def __call__(self, x, y, stderr=None, color=None, cmap='RdBu', xlabel=None, ylabel=None, title=None, label=None):
        cmap = getattr(cm, cmap)
        if color is None:
            color = 'b'
        elif type(color) != str:
            color = cmap(color)
        
        fig, ax = self._get_figure_axis(title=title, xlabel=xlabel,
                                    ylabel=ylabel)
        
        if self.ylim is not None:
            if max(y) > self.ylim[1]:
                warnings.warn('ylimit may be too small, ymax=%s' % max(y),
                        UserWarning)
            elif 2*max(y) < self.ylim[1]:
                warnings.warn('ylimit may be too large, ymax=%s' % max(y),
                        UserWarning)
            elif min(y) > self.ylim[1]:
                warnings.warn('ylimit may be too small, ymin=%s' % min(y),
                        UserWarning)
            
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
    def __call__(self, x, y_series, color_series=None, alpha=None, 
      series_labels=None, label_coords=None, cmap='RdBu', colorbar=False,
      xlabel=None, ylabel=None, title=None, filename_series=None, ui=None):
        cmap_r = getattr(cm, '%s_r' % cmap)
        cmap = getattr(cm, cmap)
        bbox = dict(facecolor='b', alpha=0.5)
        
        if color_series is not None:
            assert len(y_series) == len(color_series)
        
        if series_labels is not None:
            assert len(y_series) == len(series_labels)
            assert label_coords is not None
            label_x, label_y = label_coords
        
        fig, ax = self._get_figure_axis(title=title, xlabel=xlabel,
                                    ylabel=ylabel)
        
        if self._colorbar and colorbar:
            # probably need to set a limit on how big this will be
            ax2 = fig.add_axes([0.925, 0.1, 0.025, 0.8])
            cb = ColorbarBase(ax2, ticks=[0.0, 1.0], cmap=cmap_r,
                                        orientation='vertical')
            cb.set_ticklabels(['Low', 'High'])
            ax = fig.sca(ax) # need to make main axis the current axis again
        
        num = len(y_series) # num = window_size when plotting top genes only
        ymaxs = []
        ymins = []
        for i in ui.series(range(num), noun='Applying lines to plot'):
            if color_series is None:
                color = 'b'
            elif type(color_series[i]) != str:
                color = cmap(color_series[i])
            elif type(color_series[i]) == str:
                color = color_series[i]
            
            if series_labels is not None:
                # TODO remove hard-coded label font size
                txt = ax.text(label_x, label_y, series_labels[i],
                        bbox=bbox, color='w', fontsize=14)

            if color_series is not None:
                y = y_series[i]
            else:
                y = y_series

            ymaxs.append(max(y))
            ymins.append(min(y))

            pyplot.plot(x, y, color=color, linewidth=self.linewidth, alpha=alpha)
            
            if filename_series is not None:
                pyplot.savefig(filename_series[i])
            
            if series_labels is not None:
                ax.texts.remove(txt)

        if self.ylim is not None:
            if max(ymaxs) > self.ylim[1]:
                warnings.warn('ylimit [%s] may be too small, ymax=%s' % \
                (self.ylim[1], max(ymaxs)), UserWarning)
            elif 2*max(ymaxs) < self.ylim[1]:
                warnings.warn('ylimit may be too large, ymax=%s' % max(ymaxs),
                        UserWarning)
            elif min(ymins) > self.ylim[1]:
                warnings.warn('ylimit may be too small, ymin=%s' % min(ymins),
                        UserWarning)
        
    
