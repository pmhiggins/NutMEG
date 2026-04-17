import math
from .maintenance_model import MaintenanceModel

class Lever2015(MaintenanceModel):
    """
    Calculates cell-specific maintenance power according to Lever et al. 2015,
    specifically protein racemization. Intended to be a useful lower limit
    temperature-dependent maintenace power.

    Attributes
    ----------
    requires : dict or Nonetype
        Dictionary of additional required host properties to run
        this GrowthModel. Keys are property identifiers, and values are the
        object in a NutMEG.core class to look in.
    cutoff_pc : float
        The % racemization of AA in proteins at which it must be replaced.
    """

    def __init__(self, cutoff_pc=10):
        """
        Parameters
        ----------
        cutoff_pc : float
            The % racemization of AA in proteins at which it must be replaced.
        """
        super().__init__()
        self.cutoff_pc = cutoff_pc


    def compute(self, host, locale):
        """
        Calculate and return the ower cost due to temperature following
        Lever et al., 2015, for protein racemization.

        Parameters
        ----------
        host : base_organism
            Host organism. Must have a correct dry_mass attribute for this
            calculation to be accurate.
        locale : reactor
            Host chemical reactor. Must have the correct temperature for
            this calculation to be accurate.
        """
        if not host.E_synth:
            raise ValueError("Unable to calculate Lever 2015 maintenance power as host's E_synth is not defined")


        k_yr = 0.00012*math.exp(0.10174*(locale.env.T-273.15))
        k_s = k_yr/(365*24*3600)
        return (100*host.E_synth*k_s)/self.cutoff_pc
