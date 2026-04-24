from .forcing_factor import ForcingFactor

class Monod(ForcingFactor):
    """
    Calucates Monod-style forcing factor for microbial
    kinetic models.

    Attributes
    ----------
    substrate : str
        Name of chemical species which is rate-limiting.
    K_s : float
        Half-saturation constant with respect to substrate.
    conctype : str
        Concentration identifier of ``reagent`` to use. Can be either 'conc',
        'molal', 'activity'.
    requires : dict or Nonetype
        Dictionary of additional required host properties to run
        this GrowthModel. Keys are property identifiers, and values are the
        object in a NutMEG.core class to look in.
    """

    def __init__(self, substrate, K_s, conctype='molality'):
        """
        Extends ForcingFactor.__init__()

        Parameters
        ----------
        substrate : str
            Name of chemical species which is rate-limiting.
        K_s : float
            Half-saturation constant with respect to substrate.
        conctype : str, optional
            Concentration identifier of ``reagent`` to use. Can be either 'conc',
            'molal', 'activity'.
        """
        super().__init__()
        self.substrate = substrate
        self.K_s = K_s
        self.conctype = conctype

    def compute(self, host, locale):
        """
        Calculate and return the Monod forcing factor.

        Parameters
        ----------
        host : base_organism
            Host organism. Can be passed as None for this ForcingFactor.
        locale : reactor
            Host chemical reactor. Must have a correct concentration of
            ``substrate``.
        """
        S = None
        if self.conctype == 'molality':
            S = locale.composition[self.substrate].molality
        elif self.conctype == 'molarity':
            S = locale.composition[self.substrate].molarity
        elif self.conctype == 'activity':
            S = locale.composition[self.substrate].activity
        else:
            raise ValueError('Unknown conctype passed to Monod model')
        return S / (S + self.K_s)
