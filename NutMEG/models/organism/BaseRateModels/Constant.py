from .BaseRateModel import BaseRateModel

class Constant(BaseRateModel):
    """
    A fixed maximum rate.

    Attributes
    ----------
    requires : list
        List of additional required host properties to calucate
        this ForcingFactor.
    """

    def __init__(self, host, locale, val):
        """
        extends BaseRateModel.__init__()

        Parameters
        ----------
        host : ``base_organism'' like
            Host organism. Can be passed as None for this calculation.
        val : float
            Value of maximum rate
        """

        super().__init__(host, locale)
        self.val = val

    def compute(self, host, locale):
        """ Calculate and return the base rate for this process."""
        return self.val
