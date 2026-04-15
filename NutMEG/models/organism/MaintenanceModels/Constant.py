from .MaintenanceModel import MaintenanceModel

class Constant(MaintenanceModel):
    """
    A fixed maintneance power

    Attributes
    ----------
    requires : list
        List of additional required host properties to calucate
        this Maintenance Power.
    """

    def __init__(self, host, locale, val):
        """
        extends MaintenancePower.__init__()

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
        """ return the maintnenace power."""
        return self.val
