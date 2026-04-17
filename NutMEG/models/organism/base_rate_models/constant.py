from .base_rate_model import BaseRateModel

class Constant(BaseRateModel):
    """
    A fixed maximum rate.

    Attributes
    ----------
    requires : dict or Nonetype
        Dictionary of additional required host properties to run
        this GrowthModel. Keys are property identifiers, and values are the
        object in a NutMEG.core class to look in.
    """

    def __init__(self, val):
        """
        extends BaseRateModel.__init__()

        Parameters
        ----------
        val : float
            Value of maximum rate
        """

        super().__init__()
        self.val = val

    def compute(self, host=None, locale=None):
        """ Calculate and return the base rate for this process."""
        return self.val
