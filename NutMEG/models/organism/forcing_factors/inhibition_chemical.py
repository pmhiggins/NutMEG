from .forcing_factor import ForcingFactor

class InhibitionChemical(ForcingFactor):
    r"""
    A generic forcing factor characterizing kinetic inhibition owing to an
    inhibitory chemical species:

    .. math::

        F = \frac{1}{1 + (m/K_i)}

    where m is the concentration of the inhibitor and :math:`K_i` is its inhibition
    constant, which corresponds to the concentration at which optimum rate is
    halved.

    Attributes
    ----------
    inhibitor : str
        Name of inhibitary chemical species which is rate-limiting.
    K_i : float
        Inhibition constant with respect to inhibitor.
    conctype : str
        Concentration identifier of ``reagent`` to use. Can be either 'conc',
        'molal', 'activity'.
    requires : dict or Nonetype
        Dictionary of additional required host properties to run
        this GrowthModel. Keys are property identifiers, and values are the
        object in a NutMEG.core class to look in.
    """

    def __init__(self, inhibitor, K_i, conctype='molality',):
        """
        Extends ForcingFactor.__init__()

        Parameters
        ----------
        inhibitor : str
            Name of inhibitary chemical species which is rate-limiting.
        K_i : float
            Inhibition constant with respect to inhibitor.
        conctype : str
            Concentration identifier of ``reagent`` to use. Can be either 'conc',
            'molal', 'activity'.
        """
        super().__init__()
        self.inhibitor = inhibitor
        self.K_i = K_i
        self.conctype = conctype

    def compute(self, host, locale):
        """
        Calculate and return the inhibition forcing factor.

        Parameters
        ----------
        host : base_organism
            Host organism. Can be passed as None for this ForcingFactor.
        locale : reactor
            Host chemical reactor. Must have a correct concentration of
            ``inhibitor``.
        """
        I = locale.composition[self.inhibitor].get_amount(self.conctype)
        return 1 / (1 + (I/self.K_i))
