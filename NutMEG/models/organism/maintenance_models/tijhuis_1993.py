import math
from .maintenance_model import MaintenanceModel

class Tijhuis1993(MaintenanceModel):
    """
    Calculates cell-specific maintenance power according to Tijhuis1993

    Attributes
    ----------
    requires : dict or Nonetype
        Dictionary of additional required host properties to run
        this GrowthModel. Keys are property identifiers, and values are the
        object in a NutMEG.core class to look in.
    fit : str
        Which type of organism fit from Tijhuis 1993 to use.
    """

    def __init__(self, fit='average'):
        """
        Parameters
        ----------
        fit : str
            Which fit to use: 'average' will use the fit for an average cell,
            'anerobe' will use the fit for anaerobes, and 'aerobe' will use the
            fit for aerobes. See the Tijhuis 1993 article for details.
        """
        super().__init__()
        self.fit = fit


    def compute(self, host, locale):
        """
        Calculate and return the ower cost due to temperature following
        Tijhuis 1993

        host : base_organism
            Host organism. Must have a correct dry_mass attribute for this
            calculation to be accurate.
        locale : reactor
            Host chemical reactor. Must have the correct temperature for
            this calculation to be accurate.
        """
        if not host.dry_mass:
            raise ValueError("Unable to calculate Tijhuis 1993 maintenance power as host's dry_mass is not defined")

        drym_ME = None
        if self.fit == 'average':
            drym_ME = 4.5*math.exp((locale.env.T**(-1)-298**(-1))*((-6.94*10000)/8.31))
        elif self.fit == 'aerobe':
            5.7*math.exp((locale.env.T**(-1)-298**(-1))*((-6.94*10000)/8.31))
        elif self.fit == 'anaerobe':
            drym_ME = 3.3*math.exp((locale.env.T**(-1)-298**(-1))*((-6.94*10000)/8.31))
        else:
            return ValueError('Unknown Tijhuis et al 1993 fit')

        return (1000/3600)*(host.dry_mass/0.026)*drym_ME
