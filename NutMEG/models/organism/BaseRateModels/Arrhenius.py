from BaseRateModel import BaseRateModel

class FirstOrderChemical:
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
    E_a : float
        Activation energy
    """

    def __init__(self, host, A, Ea):
        """
        extends BaseRateModel.__init__()

        TODO: add ability to pass alternative k values, like k_RTP.

        Parameters
        ----------
        host : ``base_organism'' like
            Host organism. Can be passed as None for this initialisation.
        A : float
            Pre-exponential Arrhenius factor
        E_a : float
            Activation energy
        """

        super().__init__(host)
        self.A = A
        self.E_a = E_a


    def compute(self, host):
        """
        Calculate and return the base rate for this process, following an
        Arrhenius law. This is the default rate calculation in nutmeg.reaction.
        """

        host.net_pathway.frequency_factor = self.A
        host.net_pathway.molar_activation_E = self.E_a
        host.net_pathway.calculate_rate()
        _r = host.net_pathway.rate_const_env
        return _r
