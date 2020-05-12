import NutMEG.plotter.plotter_setup as nf
import matplotlib.pyplot as plt

class growthparams(nf.nutfig):

    def __init__(self, fig=None):
        super().__init__(fig=fig, figsize=(7,7))

    def linearplot(self, ax=None, xvals=[[]], yvals=[[]], labels=[], colors=[],
    cs='default', show=False):
        """plot the list of list values on a linear plot on the axis given.
        If no axis is passed, uses the whole figure."""
        if not len(xvals) == len(yvals):
            raise ValueError('labels, xvals and yvals do not match')

        if colors == []:
            # empty color list, either use a specific style, or the default.
            colors, lines = growthparams.colorpicker(len(xvals), cs)
        if labels == []:
            #no labels
            for x in xvals:
                labels.append('')

        if ax== None:
            # no axes have been passed, so put it in the first one
            ax = self.fig.add_subplot(111)
        for x, y, l, c, ls in zip(xvals, yvals, labels, colors, lines):
            self.addlinearplot(ax, x, y, label=l, color=c, ls=ls)
        # plt.legend()
        if show:
            plt.show()
        else:
            return ax

    def addlinearplot(self, ax, xrange, yrange, label='', color=None, ls='-'):
        """Add one linear plot to the axis ax."""
        ax.plot(xrange, yrange, label=label, c=color, linestyle=ls)


    def colorpicker(length, cs='default'):
        if cs=='default':
            # color the plot using matplotlib's standard color order
            return [None]*length, ['-']*length

        elif len(cs) ==1:
            # all one color, dash any secondary plots
            lslst=[]
            clst=[]
            for i, a in enumerate(range(length)):
                if i==0:
                    lslst.append('-')
                else:
                    lslst.append('--')
                clst.append(cs[0])
            return clst, lslst

        else:
            # request is not counted for yet, use matplotlib's standard order
            return [None]*length, ['-']*length
