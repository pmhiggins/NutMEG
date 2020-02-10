import numpy as np

import NutMEG.util.NutMEGparams as nmp
from NutMEG.util.loggersetup import loggersetup as logset
logger = logset.get_logger(__name__, filelevel=nmp.filelevel, printlevel=nmp.printlevel)


class horde_output:

    params={}

    def __init__(self, hosthorde):
        self.hosthorde = hosthorde
        params = self.refreshparams()

    def buildbaseparams(self):
        n = self.hosthorde.OrgID
        self.params = {'no_alive_'+n: np.array([]),
          'EnergyAvailable_'+n: np.array([]),
          'EnergyConservable_'+n: np.array([]),
          'CatabolicRate_'+n: np.array([]),
          'MaintenanceFrac_'+n: np.array([]),
          'Volume_'+n: np.array([]),
          'PowerSupply_'+n: np.array([]),
          'GrowthPower_'+n: np.array([]),
          'GrowthRate_'+n: np.array([]),
          'CHNOPSUptakes_'+n:np.array([])}


    def refreshparams(self):
        """reset the parameter dictionary to have empty numpy arrays.

        Inherit and adjust this method for saved organisms which have
        specific parameters you want to monitor."""
        self.buildbaseparams()

    def appendvals(self, startnum, dt):
        n = self.hosthorde.OrgID
        hc = self.hosthorde

        self.params['no_alive_'+n] = (
          np.append(self.params['no_alive_'+n],
            hc.get_population()))

        self.params['EnergyAvailable_'+n] = (
          np.append(self.params['EnergyAvailable_'+n],
            hc.respiration.net_pathway.molar_gibbs))

        self.params['EnergyConservable_'+n] = (
          np.append(self.params['EnergyConservable_'+n],
            hc.respiration.G_C))

        self.params['CatabolicRate_'+n] = (
          np.append(self.params['CatabolicRate_'+n],
            hc.respiration.rate))

        self.params['MaintenanceFrac_'+n] = (
          np.append(self.params['MaintenanceFrac_'+n],
            sum(hc.maintenance.frac_dict.values())))

        self.params['Volume_'+n] = (
          np.append(self.params['Volume_'+n],
            hc.get_volume()))

        self.params['PowerSupply_'+n] = (
          np.append(self.params['PowerSupply_'+n],
            hc.P_s))

        self.params['GrowthPower_'+n] = (
          np.append(self.params['GrowthPower_'+n],
            hc.P_growth))

        self.params['GrowthRate_'+n] = (
          np.append(self.params['GrowthRate_'+n],
            ((hc.num/startnum-1))/dt))

        # CHNOPS
        self.params['CHNOPSUptakes_'+n] = (
          np.append(self.params['CHNOPSUptakes_'+n],
            hc.CHNOPS.get_uptake(numcells=hc.num)))


        # throttled_E_growth?

        # Limiter?

        # Throttling?

    def paramtypelst(self):
        n = self.hosthorde.OrgID
        return [('no_alive_'+n, 'REAL'),
          ('EnergyAvailable_'+n, 'REAL'),
          ('EnergyConservable_'+n, 'REAL'),
          ('CatabolicRate_'+n, 'REAL'),
          ('MaintenanceFrac_'+n, 'REAL'),
          ('Volume_'+n, 'REAL'),
          ('PowerSupply_'+n, 'REAL'),
          ('GrowthPower_'+n, 'REAL'),
          ('GrowthRate_'+n, 'REAL'),
          ('CHNOPSUptakes_'+n, 'TEXT')]
