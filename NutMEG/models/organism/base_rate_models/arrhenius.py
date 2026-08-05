from .base_rate_model import BaseRateModel
from NutMEG.utils.math import math

class Arrhenius(BaseRateModel):
    """
    A maximum rate calculated using a first-order chemical calculation
    with the net metabolic pathway.

    Attributes
    ----------
    requires : list
        List of additional required host properties to calucate
        this ForcingFactor.
    A : float
        Pre-exponential Arrhenius factor
    Ea : float
        Molar activation energy. Unit J/K mol.
    requires : dict or Nonetype
        Dictionary of additional required host properties to run
        this GrowthModel. Keys are property identifiers, and values are the
        object in a NutMEG.core class to look in.
    """

    def __init__(self, A, Ea, **kwargs):
        """
        extends BaseRateModel.__init__()

        TODO: add ability to have E_a as a temperature-dependent function.

        Parameters
        ----------
        A : float
            Pre-exponential Arrhenius factor
        Ea : float
            Molar activation energy. Unit J/K mol.
        """

        super().__init__()
        self.A = A
        self.Ea = Ea


    def compute(self, host, locale):
        """
        Calculate and return the base rate for this process, following an
        Arrhenius law.

        Parameters
        ----------
        host : base_organism or None
            Host organism. Can be passed as None for this BaseRateModel.
        locale : reactor
            Host chemical reactor. Must have the correct temperature for
            this calculation to be accurate.
        """

        return self.A * math.exp(-self.Ea/(9.81*locale.T))
