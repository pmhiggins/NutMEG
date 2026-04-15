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

    def __init__(self, host, val):
        """
        extends BaseRateModel.__init__()

        Parameters
        ----------
        host : ``base_organism'' like
            Host organism. Can be passed as None for this calculation.
        val : float
            Value of maximum rate
        """

        super().__init__()
        self.val = val

    def compute(self, host):
        """ Calculate and return the base rate for this process."""
        return self.val
