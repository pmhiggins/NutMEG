from .ForcingFactor import ForcingFactor

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
    """

    def __init__(self, host, locale, substrate, K_s, conctype='molal'):
        """
        Extends ForcingFactor.__init__()

        Parameters
        ----------
        host : ``base_organism'' like
            Host organism. Not strictly required for this ForcingFactor
            instance, so may be passed as None
        substrate : str
            Name of chemical species which is rate-limiting.
        K_s : float
            Half-saturation constant with respect to substrate.
        """
        super().__init__(host, locale)
        self.substrate = substrate
        self.K_s = K_s
        self.conctype = conctype

    def compute(self, host, locale):
        S = None
        if self.conctype == 'molal':
            S = locale.composition[self.substrate].molal
        elif self.conctype == 'conc':
            S = locale.composition[self.substrate].conc
        elif self.conctype == 'activity':
            S = locale.composition[self.substrate].activity
        else:
            raise ValueError('Unknown conctype passed to Monod model')
        return S / (S + self.K_s)
