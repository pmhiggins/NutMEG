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
    A : float
        Pre-exponential Arrhenius factor
    E_a : float
        Activation energy
    requires : dict or Nonetype
        Dictionary of additional required host properties to run
        this GrowthModel. Keys are property identifiers, and values are the
        object in a NutMEG.core class to look in.
    """

    def __init__(self, A, Ea):
        """
        extends BaseRateModel.__init__()

        TODO: add ability to pass alternative k values, like k_RTP.

        Parameters
        ----------
        A : float
            Pre-exponential Arrhenius factor
        E_a : float
            Activation energy
        """

        super().__init__()
        self.A = A
        self.E_a = E_a


    def compute(self, host, locale):
        """
        Calculate and return the base rate for this process, following an
        Arrhenius law. This is the default rate calculation in NutMEG.reaction.

        The passed host.metabolism.net_pathway will have its rate constants,
        activation energies and pre-exponential factors updated based on the
        attributes of this object.

        Parameters
        ----------
        host : base_organism
            Host organism.
        locale : reactor
            Host chemical reactor. Must have the correct temperature for
            this calculation to be accurate.
        """

        host.metabolism.net_pathway.frequency_factor = self.A
        host.metabolism.net_pathway.molar_activation_E = self.E_a
        host.metabolism.net_pathway.calculate_rate()
        _r = host.metabolism.net_pathway.rate_const_env
        return _r
