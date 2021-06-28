import numpy as np

import NutMEG.util.NutMEGparams as nmp
from NutMEG.util.loggersetup import loggersetup as logset
logger = logset.get_logger(__name__, filelevel=nmp.filelevel, printlevel=nmp.printlevel)


class horde_output:
    """helper class for managing the output of a horde, either to terminal
    or to a database.

    Attributes
    ----------
    hosthorde : horde
        The `host` for this output generator, Where we get our numbers from.
    params : dict
        Dictionary of parameters to output.
    """
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
          'CHNOPSUptakes_'+n:np.array([]),
          'tot_no_cells_'+n:np.array([]),
          'Volume_tot_'+n:np.array([])}


    def refreshparams(self):
        """reset the parameter dictionary to have empty numpy arrays.

        Inherit and adjust this method for saved organisms which have
        specific parameters you want to monitor."""
        self.buildbaseparams()

    def appendvals(self, startnum, dt):
        """Add parameters at the present time to params.

        Parameters
        ----------
        startnum : int
            `Previous` number of organisms, assuming this is output after a
            time step. Used to calculate thr gorwth rate
        dt : float
            Amount of time that has passed since the previous step.
        """
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

        self.params['tot_no_cells_'+n] = (
          np.append(self.params['tot_no_cells_'+n],
            hc.get_population(inactive=True)))

        self.params['Volume_tot_'+n] = (
          np.append(self.params['Volume_tot_'+n],
            hc.get_volume(inactive=True)))


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
          ('CHNOPSUptakes_'+n, 'TEXT'),
          ('tot_no_cells_'+n, 'REAL'),
          ('Volume_tot_'+n, 'REAL')]
