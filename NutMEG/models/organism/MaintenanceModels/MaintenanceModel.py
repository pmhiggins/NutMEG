

class MaintenanceModel:
    """
    Superclass for all maintenance model calculations.

    Attributes
    ----------
    requires : dict or Nonetype
        Dictionary of additional required host properties to run
        this GrowthModel. Keys are property identifiers, and values are the
        object in a NutMEG.core class to look in.
    """

    def __init__(self):

        self.requires = None

    def compute(self, host, locale):
        """ Calculate and return the ForcingFactor for this process."""
        raise NotImplementedError

    def outputs(self):
        """
        Return a dict of the key outputs for host properties this calculation
        generated.
        """
        return {}
