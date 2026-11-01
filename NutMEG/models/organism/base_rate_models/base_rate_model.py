


class BaseRateModel:
    """
    Superclass for all base rate calculations (i.e., methods to calculate
    maximum rate constants or growth rates).

    Attributes
    ----------
    requires : dict or Nonetype
        Dictionary of additional required host properties to run
        this GrowthModel. Keys are property identifiers, and values are the
        object in a NutMEG.core class to look in.
    """

    def __init__(self, **kwargs):

        self.requires = None

    def compute(self, host, locale):
        """ Calculate and return the base rate for this process."""
        raise NotImplementedError("A max rate likely hasn't been defined.")

    def outputs(self):
        """
        Return a dict of the key outputs for host properties this calculation generated.
        """
        return {}
