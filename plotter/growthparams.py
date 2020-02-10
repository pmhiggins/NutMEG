import NutMEG.plotter.plotter_setup as nf
import matplotlib.pyplot as plt

class growthparams(nf.nutfig):

    def __init__(self, fig=None):
        super().__init__(fig=fig, figsize=(7,7))

    def linearplot(self, ax=None, xvals=[[]], yvals=[[]], labels=[], colors=[], show=False, ls='-'):
        """plot the list of list values on a linear plot on the axis given.
        If no axis is passed, uses the whole figure."""
        if not len(xvals) == len(yvals) == len(labels):
            raise ValueError('labels, xvals and yvals do not match')
        if ax== None:
            # no axes have been passed, so put it in the first one
            ax = self.fig.add_subplot(111)
        for x, y, l, c in zip(xvals, yvals, labels, colors):
            self.addlinearplot(ax, x, y, label=l, color=c, ls=ls)
        # plt.legend()
        if show:
            plt.show()
        else:
            return ax

    def addlinearplot(self, ax, xrange, yrange, label='', color=None, ls='-'):
        """Add one linear plot to the axis ax."""
        ax.plot(xrange, yrange, label=label, c=color, linestyle=ls)
