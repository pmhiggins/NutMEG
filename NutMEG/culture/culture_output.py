import numpy as np

import NutMEG.util.NutMEGparams as nmp
from NutMEG.util.loggersetup import loggersetup as logset
logger = logset.get_logger(__name__, filelevel=nmp.filelevel, printlevel=nmp.printlevel)

class culture_output:

    params={}

    def __init__(self, hostculture):
        self.hostculture = hostculture
        params = self.refreshparams()

    def buildbaseparams(self):
        self.params = {'totBM_cells': np.array([]),
          'totBM_vol': np.array([]),
          'totBM_kg': np.array([])}

    def refreshparams(self):
        """reset the parameter dictionary to have empty numpy arrays.
        """
        self.buildbaseparams()


    def appendvals(self):
        hc = self.hostculture

        self.params['totBM_cells'] = (
          np.append(self.params['totBM_cells'],
            hc.get_population()))

        self.params['totBM_vol'] = (
          np.append(self.params['totBM_vol'],
            hc.get_total_volume()))

        self.params['totBM_kg'] = (
          np.append(self.params['totBM_kg'],
            hc.get_total_mass()))

    def paramtypelst(self):
        return [('totBM_cells', 'REAL'),
          ('totBM_vol', 'REAL'),
          ('totBM_kg', 'REAL')]
