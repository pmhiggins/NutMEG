from .BaseRateModel import BaseRateModel

class FirstOrderChemical:
    """
    A maximum rate calculated using a first-order chemical calculation
    with the net metabolic pathway.

    Attributes
    ----------
    requires : list
        List of additional required host properties to calucate
        this ForcingFactor.
    k_env : float
        First-order rate constant in current environment condition
    """

    def __init__(self, host, k_env):
        """
        extends BaseRateModel.__init__()

        TODO: add ability to pass alternative k values, like k_RTP.

        Parameters
        ----------
        host : ``base_organism'' like
            Host organism. Can be passed as None for this calculation.
        k_env : float
            First-order rate constant in current environment condition
        """

        super().__init__(host)
        self.k_env = k_env


    def compute(self, host):
        """
        Calculate and return the base rate for this process, following a
        generic first order rate law.
        """
        conc_multiplier = 1.0
        for r, mr in host.metabolism.net_pathway.reactants.items():
            conc_multiplier = conc_multiplier*(r.activity**mr)

        return self.k_env * conc_multiplier
