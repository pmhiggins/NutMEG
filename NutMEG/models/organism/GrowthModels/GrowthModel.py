

class GrowthModel:
    """
    Superclass for all growth rate calculations.

    Attributes
    ----------
    requires : list
        List of additional required host properties to run
        this GrowthModel.
    """

    def __init__(self, host, locale):
        """
        Parameters
        ----------
        host : ``base_organism'' like
            Host organism. Some GrowthModels will need this to initialise and
            some won't. It is best to assume they will (else they may throw an error)
        locale : ``reactor'' like
            Host chemical reactor. Some GrowthModels will need this to initialise and
            some won't. It is best to assume they will (else they may throw an error)
        """

        self.requires = None

    def compute(self, host, locale):
        """ Calculate and return the growth rate according to this model."""
        raise NotImplementedError

    def outputs(self):
        """
        Return a dict of the key outputs for host properties this calculation generated.
        """
        return {}

    def find_subparams(self):
        """

        """
        subparams = {}
        for n, loc in self.requires.items():
            if type(loc) == list:
                for h in loc:
                    if n in h.outputs():
                        subparams[n] = h.outputs()[n]
            else:
                if n in h.outputs:
                    subparams[n] = h.outputs[n]
        return subparams
