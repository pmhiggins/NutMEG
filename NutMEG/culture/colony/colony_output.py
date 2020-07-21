import numpy as np

class colony_output:
    """helper class for managing the output of a colony, either to terminal
    or to a database."""

    params={}

    def __init__(self, hostcolony):
        self.hostcolony = hostcolony
        params = self.refreshparams()

    def buildbaseparams(self, first=False):
        n = self.hostcolony.OrgID
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
        n = self.hostcolony.OrgID
        hc = self.hostcolony

        self.params['no_alive_'+n] = (
          np.append(self.params['no_alive_'+n],
            hc.get_population()))

        # not ideal, average?
        self.params['EnergyAvailable_'+n] = (
          np.append(self.params['EnergyAvailable_'+n],
            hc.collection[0].respiration.net_pathway.molar_gibbs))

        # not ideal, average?
        self.params['EnergyConservable_'+n] = (
          np.append(self.params['EnergyConservable_'+n],
            hc.collection[0].respiration.G_C))

        # not ideal, average?
        self.params['CatabolicRate_'+n] = (
          np.append(self.params['CatabolicRate_'+n],
            hc.collection[0].respiration.rate))

        self.params['MaintenanceFrac_'+n] = (
          np.append(self.params['MaintenanceFrac_'+n],
            hc.get_average_maintenance_fraction()))

        self.params['Volume_'+n] = (
          np.append(self.params['Volume_'+n],
            hc.get_volume()))

        # not ideal. Average?
        self.params['PowerSupply_'+n] = (
          np.append(self.params['PowerSupply_'+n],
            hc.base_Ps))

        # not ideal. Average?
        self.params['GrowthPower_'+n] = (
          np.append(self.params['GrowthPower_'+n],
            hc.collection[0].P_growth))

        self.params['GrowthRate_'+n] = (
          np.append(self.params['GrowthRate_'+n],
            ((len(hc.collection)/startnum)-1)/dt))

        # CHNOPS
        self.params['CHNOPSUptakes_'+n] = (
          np.append(self.params['CHNOPSUptakes_'+n],
            hc.getCHNOPSut()))


        # throttled_E_growth?

        # Limiter?

        # Throttling?


    def paramtypelst(self):
        n = self.hostcolony.OrgID
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
