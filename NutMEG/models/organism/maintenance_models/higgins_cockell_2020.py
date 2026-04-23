import math
from .maintenance_model import MaintenanceModel

class HigginsCockell2020(MaintenanceModel):
    """
    Calculates cell-specific maintenance power according to Higgins and Cockell
    2020. Specifically a fit for hydrogenotrophic methanogens at varying
    temperatures.

    Attributes
    ----------
    requires : dict or Nonetype
        Dictionary of additional required host properties to run
        this GrowthModel. Keys are property identifiers, and values are the
        object in a NutMEG.core class to look in.
    """

    # polyfitting parameters from Higgins and Cockell 2020 at different ATP yields.
    TOM_1 = [ 3.75274464e-10, -6.28614763e-07,  4.20934430e-04,
      -1.40873711e-01, 2.36260582e+01, -1.60782148e+03]
    TOM_05 = [ 5.47545848e-10, -9.36824600e-07,  6.41040964e-04,
      -2.19319054e-01, 3.75813051e+01, -2.59970377e+03]
    TOM_15 = [ 1.40807504e-10, -2.37559221e-07,  1.60324386e-04,
      -5.41193986e-02, 9.19861608e+00, -6.48600497e+02]

    def __init__(self, n_ATP=None):
        """
        Parameters
        ----------
        n_ATP : float
            ATP yield of the overall metabolism.
            Below or above 0.5 or 1.5, the fit is assigned to use them respectively.
            If None, then use the value fitted at n_ATP = 1. Default None.
        """
        super().__init__()
        self.n_ATP = n_ATP
        if n_ATP == None:
            self.n_ATP = 1.

    def compute(self, host, locale):
        """
        Calculate and return the power cost due to temperature following
        Higgins and Cockell 2020

        Parameters
        ----------
        host : base_organism
            Host organism. Must have a correct dry_mass attribute for this
            calculation to be accurate.
        locale : reactor
            Host chemical reactor. Must have the correct temperature for
            this calculation to be accurate.
        """
        if not host.dry_mass:
            raise ValueError("Unable to calculate HigginsCockell maintenance power as host's dry_mass is not defined")

        HCpolyfit = None

        if self.n_ATP >= 1.5:
            HCpolyfit = self.TOM_15
        elif self.n_ATP <= 0.5:
            HCpolyfit = self.TOM_05
        else:
            HCpolyfit = self.TOM_1

        MP = 10**sum([x*(locale.T**(5-i)) for i, x in enumerate(HCpolyfit)])
        return MP * host.dry_mass / (300*3.4412868852915668e-18 )
