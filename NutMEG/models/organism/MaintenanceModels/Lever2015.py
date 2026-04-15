import math
from .MaintenanceModel import MaintenanceModel

class Lever2015(MaintenanceModel):
    """
    Calculates cell-specific maintenance power according to Lever et al. 2015,
    specifically protein racemization. Intended to be a useful lower limit
    temperature-dependent maintenace power.

    Attributes
    ----------
    requires : list
        List of additional required host properties to run
        this MaintenanceModel.
    cutoff_pc : float
        The % racemization of AA in proteins at which it must be replaced.
    """

    def __init__(self, host, locale, cutoff_pc=10):
        """
        Parameters
        ----------
        host : ``base_organism'' like
            Host organism. Some MaintenanceModels will need this to initialise and
            some won't. It is best to assume they will (else they may throw an error)
        locale : ``reactor'' like
            Host chemical reactor. Some Forcing Factors will need this to initialise and
            some won't. It is best to assume they will (else they may throw an error)
        cutoff_pc : float
            The % racemization of AA in proteins at which it must be replaced.
        """
        super().__init__(host, locale)
        self.cutoff_pc = cutoff_pc


    def compute(self, host, locale):
        """
        Calculate and return the ower cost due to temperature following
        Lever et al., 2015, for protein racemization.
        """
        if not host.E_synth:
            raise ValueError("Unable to calculate Lever 2015 maintenance power as host's E_synth is not defined")

        """Power cost due to temperature according to Lever et al 2015
        """

        k_yr = 0.00012*math.exp(0.10174*(locale.env.T-273.15))
        k_s = k_yr/(365*24*3600)
        return (100*host.E_synth*k_s)/self.cutoff_pc


    def outputs(self):
        """
        Return a dict of the key outputs for host properties this calculation generated.
        """
        return {}
