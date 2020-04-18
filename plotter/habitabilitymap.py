import NutMEG.plotter.plotter_setup as nf
import matplotlib.pyplot as plt
import sqlite3
import os
import math
import matplotlib as mpl
from matplotlib import gridspec
import numpy as np

std_dbpath=os.path.join(os.path.dirname(__file__), '../testDB')

class habitabilitymap(nf.nutfig):

    def __init__(self, fig=None):
        super().__init__(fig=fig, figsize=(7,7))

    def add_map_plot(self, ax=None, xvals=[], yvals=[], zvals=[], labels=['','',''], show=False, save='', colorbar=True):
        """pass 1D lists of the x y and z coordinates of the proposed
        map. This map will then be placed onto ax, or if you have no
        ax preference will take up the whole figure.

        If you pass colorbar as true, one will be generated in a
        seperate matplotlib figure.

        if you pass save, the Figure and/or colorbar will be saved to the
        filename given, colorbar will have the prefix cb_. Please do
        include the file extension.
        """
        if not len(xvals) == len(yvals) == len(zvals):
            raise ValueError('xvals, yvals and zvals do not match')
        if ax== None:
            # no axes have been passed, so put it in the first one
            ax = self.fig.add_subplot(111)

        # normalise x and y to get a standard plot (may not be neccessary)
        xvalsn, xma, xmi = nf.nutfig.normalise(xvals)
        yvalsn, yma, ymi = nf.nutfig.normalise(yvals)

        ax.set_xlabel(labels[0])
        ax.set_ylabel(labels[1])

        ax.tricontourf(xvalsn, yvalsn, np.array(zvals), 100, vmin=min(zvals), vmax=max(zvals), cmap='gist_earth_r')

        xa=ax.get_xticks().tolist()
        xa = [round(nf.nutfig.denormalise(aa, xma, xmi),3) for aa in xa]
        ax.set_xticklabels(xa)
        ya=ax.get_yticks().tolist()
        ya = [round(nf.nutfig.denormalise(aa, yma, ymi),3) for aa in ya]
        ax.set_yticklabels(ya)
        self.mapcolorbarh(zvals, labels[2], save=save)
        if not save == '':
            plt.savefig(save)
        if show:
            plt.show()

    def mapcolorbarh(self, zvals, zlabel, save=''):
        """Produce a horzontal colorbar on a new canvas with the zvals and
        label passed."""
        cbfig = plt.figure(figsize=(0.0393701*100, 0.0393701*40))
        gs = gridspec.GridSpec(1, 1, left=0.05, right=0.95, top=0.9, bottom=0.55)


        cbaxes = cbfig.add_subplot(gs[0])
        # cbaxes.tick_params(axis='both', which='major', labelsize=16)
        norm=mpl.colors.Normalize(vmin=min(zvals), vmax=max(zvals))
        cb2 = mpl.colorbar.ColorbarBase(cbaxes,
                            cmap=plt.get_cmap('gist_earth_r'),
                            norm=norm,
                            orientation='horizontal')
        cb2.set_label(zlabel)#, fontsize=16)

        a=cbaxes.get_xticks().tolist()
        a = ['{:0.1e}'.format(aa) for aa in a]

        cbaxes.set_xticklabels(a)
        if not save == '':
            plt.savefig('cb_'+save)

    def mapdata_orgloc(self, xIDs, yIDs, zname):
        """ for the OrgIDs and LocIDs passed, this finds the value of the zname parameter in the database corresponding to each.
        """
        db=sqlite3.connect(std_dbpath)
        # db.row_factory = lambda cursor, row: row[0]
        cursor = db.cursor()
        zparams=[]
        try:
            for x,y in zip(xIDs, yIDs):
                zparams.append(nf.nutfig.extract_param_db_OrgLoc(x, y, zname))
            return zparams
        except sqlite3.Error as e:
            raise e
        finally:
            db.close()
