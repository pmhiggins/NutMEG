

class MaintenanceModel:
    """
    Superclass for all maintenance model calculations.

    Attributes
    ----------
    requires : list
        List of additional required host properties to run
        this MaintenanceModel.
    """

    def __init__(self, host, locale):
        """
        Parameters
        ----------
        host : ``base_organism'' like
            Host organism. Some MaintenanceModels will need this to initialise and
            some won't. It is best to assume they will (else they may throw an error)
        locale : ``reactor'' like
            Host chemical reactor. Some MaintenanceModels will need this to initialise and
            some won't. It is best to assume they will (else they may throw an error)
        """

        self.requires = None

    def compute(self, host, locale):
        """ Calculate and return the ForcingFactor for this process."""
        return NotImplementedError

    def outputs(self):
        """
        Return a dict of the key outputs for host properties this calculation generated.
        """
        return {}
