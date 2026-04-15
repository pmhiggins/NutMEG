


class ForcingFactor:
    """
    Superclass for all forcing factor calculations.

    Attributes
    ----------
    requires : list
        List of additional required host properties to calucate
        this ForcingFactor.
    """

    def __init__(self, host, locale):
        """
        Parameters
        ----------
        host : ``base_organism'' like
            Host organism. Some Forcing Factors will need this to initialise and
            some won't. It is best to assume they will (else they may throw an error)
        """

        self.requires = None

    def compute(self, host, locale):
        """ Calculate and return the ForcingFactor for this process."""
        raise NotImplementedError

    def outputs(self):
        """
        Return a dict of the key outputs for host properties this calculation generated.
        """
        return {}
