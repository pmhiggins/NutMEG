from .base_rate_model import BaseRateModel

class FirstOrderChemical:
    """
    A maximum rate calculated using a first-order chemical calculation
    with the net metabolic pathway.

    Attributes
    ----------
    requires : dict or Nonetype
        Dictionary of additional required host properties to run
        this GrowthModel. Keys are property identifiers, and values are the
        object in a NutMEG.core class to look in.
    k_env : float
        First-order rate constant in current environment condition
    """

    def __init__(self, k_env):
        """
        extends BaseRateModel.__init__()

        TODO: add ability to pass alternative k values, like k_RTP.

        Parameters
        ----------
        k_env : float
            First-order rate constant in current environment condition
        """

        super().__init__()
        self.k_env = k_env


    def compute(self, host, locale=None):
        """
        Calculate and return the base rate for this process, following a
        generic first order rate law.

        Parameters
        ----------
        host : base_organism
            Host organism. Can be passed as None for this calculation.
        locale : reactor, optional
            Local reactor. Not currently implemented in this method.
        """
        conc_multiplier = 1.0
        for r, mr in host.metabolism.net_pathway.reactants.items():
            conc_multiplier = conc_multiplier*(r.activity**mr)

        return self.k_env * conc_multiplier
