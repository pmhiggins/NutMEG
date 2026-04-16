from .MaintenanceModel import MaintenanceModel

class Constant(MaintenanceModel):
    """
    A fixed maintneance power

    Attributes
    ----------
    requires : dict or Nonetype
        Dictionary of additional required host properties to run
        this GrowthModel. Keys are property identifiers, and values are the
        object in a NutMEG.core class to look in.
    """

    def __init__(self, val):
        """
        extends MaintenancePower.__init__()

        Parameters
        ----------
        val : float
            Value of maximum rate
        """

        super().__init__()
        self.val = val

    def compute(self, host=None, locale=None):
        """ return the maintnenace power."""
        return self.val
