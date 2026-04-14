


class BaseRateModel:
    """
    Superclass for all base rate calculations (i.e., methods to calculate
    maximum rate constants or growth rates).

    Attributes
    ----------
    requires : list
        List of additional required host properties to calucate
        this ForcingFactor.
    """

    def __init__(self, host):
        """
        Parameters
        ----------
        host : ``base_organism'' like
            Host organism. Some base rates will need this to initialise and
            some won't. It is best to assume they will (else they may throw an error)
        """

        self.requires = None

    def compute(self, host):
        """ Calculate and return the base rate for this process."""
        return NotImplementedError

    def outputs(self):
        """
        Return a dict of the key outputs for host properties this calculation generated.
        """
        return {}
