from .base_rate_model import BaseRateModel

class Q10(BaseRateModel):
    """
    A temperature-dependent maximum rate calculated using a Q10 coeffiecient and
    a reference temperature.

    .. math::

        r_T = Q^{(T-T_{ref})/10}

    It is common in biology to approximate :math:`Q = 2` to 3 and
    :math:`T_{ref}=298`. They can be adjusted to any value here. Depending on
    the specificity of the target calculation, more flexible rate-temperature laws
    may be desirable (such as a normalized biological performance function).

    Attributes
    ----------
    rate_Tref : float
        rate value at reference temperature Tref
    Q : float, optional
        Q10 coefficient, i.e. the amount by which rate_Tref increases per 10 K
        increment. Default 2.
    Tref : float, optional
        Reference temperature in K. Default 298.
    """
    def __init__(self, rate_Tref, Q=2., Tref=298.):
        """
        extends BaseRateModel.__init__()

        Parameters
        ----------
        rate_Tref : float
            rate value at reference temperature Tref
        Q : float, optional
            Q10 coefficient, i.e. the amount by which rate_Tref increases per 10 K
            increment. Default 2.
        Tref : float, optional
            Reference temperature in K. Default 298.
        """
        super().__init__()
        self.rate_Tref = rate_Tref
        self.Q = Q
        self.Tref = Tref

    def compute(self, host, locale):
        """
        Calculate and return the base rate for this process at the local temperature,
        following the stored Q10 coefficient.

        Parameters
        ----------
        host : base_organism
            Host organism. Can be passed as None for this calculation.
        locale : reactor, optional
            Local reactor.
        """
        return self.rate_Tref * (self.Q **((locale.T-self.Tref)/10.))
