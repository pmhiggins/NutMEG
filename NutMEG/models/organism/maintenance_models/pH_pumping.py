import math
from .maintenance_model import MaintenanceModel

class pH_Pumping(MaintenanceModel):
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
    pH_interior : float, optional
        pH inside the cell. Default 7.
    memb_pot : float, optional
        Membrane potential of the cell wall. Default 1e-5 V
    PermH : float, optional
        Permeability of the cell wall to Hydrogen in m^-1. Default 1e-10
    PermOH : float, optional
        Permeability of the cell wall to Hydroxide in m^-1. Default 1e-10
    """

    def __init__(self,
      membrane_potential=1e-5,
      pH_interior=7.,
      membrane_permH=1e-10,
      membrane_permOH=1e-10):
        """
        Parameters
        ----------
        pH_interior : float, optional
            pH inside the cell. Default 7.
        memb_pot : float, optional
            Membrane potential of the cell wall. Default 1e-5 V
        PermH : float, optional
            Permeability of the cell wall to Hydrogen in m^-1. Default 1e-10
        PermOH : float, optional
            Permeability of the cell wall to Hydroxide in m^-1. Default 1e-10
        """
        super().__init__()
        self.membrane_potential = membrane_potential
        self.pH_interior = pH_interior
        self.membrane_permH = membrane_permH
        self.membrane_permOH = membrane_permOH



    def compute(self, host, locale):
        """
        Calculate and return the power cost due to pH following
        the proton pumping calculation described in Higgins (2022) U. Edinburgh.

        Parameters
        ----------
        host : base_organism
            Host organism. Must have a correct surfacearea attribute for this
            calculation to be accurate.
        locale : reactor like
            Host chemical reactor. Must have the correct temperature and pH or
            [H+] for this calculation to be accurate.
        """
        if not host.surfacearea:
            raise ValueError("Unable to calculate pH_Pumping maintenance power as host's surfacearea is not defined")
        if not locale.pH:
            try:
                locale.pH = 10**-locale.composition['H+'].conc
            except:
                raise ValueError("Reactor pH unknown to pH_Pumping")

        fluxH = self._getfluxH(host, locale)
        fluxOH = self._getfluxOH(host, locale)
        E = self._getEnergyPump(host, locale)
        return (abs(fluxOH)+abs(fluxH))*E



    @staticmethod
    def GoldmanEQ(P, K, concs, SA, T, q=1):
        return P*(q)*(T/298)*(K/(1-math.exp(-q*K)))*(concs['in']-(concs['out']*math.exp(-q*K)))*SA*1e3


    def _getfluxH(self, host, locale):
        K = self.membrane_potential*96485/(8.31*locale.env.T)
        concs = {'in':10**-self.pH_interior, 'out':10**-locale.pH}
        return pH_Pumping.GoldmanEQ(self.membrane_permH, K, concs, host.surfacearea, locale.env.T, q=1)

    def _getfluxOH(self, host, locale):
        K = self.membrane_potential*96485/(8.31*locale.env.T)
        concs = {'in':10**(self.pH_interior-14),
          'out':10**(locale.pH-14)}
        return pH_Pumping.GoldmanEQ(self.membrane_permOH, K, concs, host.surfacearea, locale.env.T, q=-1)

    def _getEnergyPump(self, host, locale):
        return(8.31*locale.env.T*abs(self.pH_interior-locale.pH))
