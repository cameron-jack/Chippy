from __future__ import division

from matplotlib import pyplot, rc, cm
from matplotlib.ticker import MultipleLocator

from cogent.util.progress_display import display_wrap

class _Plottable(object):
    """base class for handling plotting"""
    def __init__(self, height, width, bgcolor, grid, pad=10, ylim=None,
                xlim=None, xtick_space=None, ytick_space=None,
                xtick_interval=None, ytick_interval=None, linewidth=2,
                xlabel_fontsize=None, ylabel_fontsize=None, vline=None,
                ioff=None):
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
    
    def _get_figure_axis(self, title=None, xlabel=None, ylabel=None):
        """returns the figure and axis ready for display"""
        if self.xlabel_fontsize:
            rc('xtick', labelsize=self.xlabel_fontsize)
        if self.ylabel_fontsize:
            rc('ytick', labelsize=self.ylabel_fontsize)
        
        fig = pyplot.figure(figsize=(self.width, self.height))
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
            pyplot.ylabel(ylabel)
        
        if xlabel:
            pyplot.xlabel(xlabel)
        
        return fig, ax

class PlottableSingle(_Plottable):
    """Plots a single line"""
    def __init__(self, *args, **kwargs):
        super(PlottableSingle, self).__init__(*args, **kwargs)
    
    def __call__(self, x, y, stderr=None, color=None, cmap='RdBu', xlabel=None, ylabel=None, title=None, filename=None):
        cmap = getattr(cm, cmap)
        if color is None:
            color = 'b'
        elif type(color) != str:
            color = cmap(color)
        
        fig, ax = self._get_figure_axis(title=title, xlabel=xlabel,
                                    ylabel=ylabel)
        pyplot.plot(x, y, color=color, linewidth=self.linewidth)
        
        if stderr is not None:
            upper = 1.96 * stderr + y
            lower = -1.96 * stderr + y
            pyplot.fill_between(x, upper, lower, alpha=0.2, color=color)
        
        if filename is not None:
            pyplot.savefig(filename)
    

class PlottableGroups(_Plottable):
    """plot groups of data on the same panel"""
    def __init__(self, *args, **kwargs):
        super(PlottableGroups, self).__init__(*args, **kwargs)
    
    @display_wrap
    def __call__(self, x, y_series, color_series=None, alpha=None, 
      series_labels=None, label_coords=None, cmap='RdBu', xlabel=None,
      ylabel=None, title=None, filename_series=None, filename=None,
      show_final=True, ui=None):
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
        num = len(y_series)
        for i in ui.series(range(num), noun='Applying lines to plot'):
            if color_series is None:
                color = 'b'
            elif type(color_series[i]) != str:
                color = cmap(color_series[i])
            
            if series_labels is not None:
                # TODO remove hard-coded label font size
                txt = ax.text(series_labels[i], label_x, label_y,
                        bbox=bbox, color='w', fontsize=14)
            
            y = y_series[i]
            pyplot.plot(x, y, color=color, linewidth=self.linewidth,
                    alpha=alpha)
            if filename_series is not None:
                pyplot.savefig(filename_series[i])
            
            if series_labels is not None:
                ax.texts.remove(txt)
        
        pyplot.ion()
        if show_final:
            pyplot.show()
        
        if filename is not None:
            pyplot.savefig(filename)
    
