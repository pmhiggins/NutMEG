import math
from .MaintenanceModel import MaintenanceModel

class pH_Pumping(MaintenanceModel):
    """
    Calculates cell-specific maintenance power according to Lever et al. 2015,
    specifically protein racemization. Intended to be a useful lower limit
    temperature-dependent maintenace power.

    Attributes
    ----------
    requires : list
        List of additional required host properties to run
        this MaintenanceModel.
    pH_interior : float, optional
        pH inside the cell. Default 7.
    memb_pot : float, optional
        Membrane potential of the cell wall. Default 1e-5 V
    PermH : float, optional
        Permeability of the cell wall to Hydrogen in m^-1. Default 1e-10
    PermOH : float, optional
        Permeability of the cell wall to Hydroxide in m^-1. Default 1e-10
    """

    def __init__(self, host, locale,
      membrane_potential=1e-5,
      pH_interior=7.,
      membrane_permH=1e-10,
      membrane_permOH=1e-10):
        """
        Parameters
        ----------
        host : ``base_organism'' like
            Host organism. Some MaintenanceModels will need this to initialise and
            some won't. It is best to assume they will (else they may throw an error)
        locale : ``reactor'' like
            Host chemical reactor. Some Forcing Factors will need this to initialise and
            some won't. It is best to assume they will (else they may throw an error)
        pH_interior : float, optional
            pH inside the cell. Default 7.
        memb_pot : float, optional
            Membrane potential of the cell wall. Default 1e-5 V
        PermH : float, optional
            Permeability of the cell wall to Hydrogen in m^-1. Default 1e-10
        PermOH : float, optional
            Permeability of the cell wall to Hydroxide in m^-1. Default 1e-10
        """
        super().__init__(host, locale)
        self.membrane_potential = membrane_potential
        self.pH_interior = pH_interior
        self.membrane_permH = membrane_permH
        self.membrane_permOH = membrane_permOH



    def compute(self, host, locale):
        """
        Calculate and return the ower cost due to temperature following
        Lever et al., 2015, for protein racemization.
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


    def outputs(self):
        """
        Return a dict of the key outputs for host properties this calculation generated.
        """
        return {}


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
