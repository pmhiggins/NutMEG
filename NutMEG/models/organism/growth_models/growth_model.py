

class GrowthModel:
    """
    Superclass for all growth rate calculations.

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
        """
        Calculate and return the growth rate according to this model.

        Parameters
        ----------
        host : base_organism
            Host organism. Some GrowthModels will need this to initialise and
            some won't. It is best to assume they will (else they may throw an error)
        locale : reactor
            Host chemical reactor. Some GrowthModels will need this to initialise and
            some won't. It is best to assume they will (else they may throw an error)
        """
        raise NotImplementedError

    def outputs(self):
        """
        Return a dict of the key outputs for host properties this calculation
        generated. Its format is an identifier key, followed by the attribute.
        """
        return {}

    def find_subparams(self, alt_requires=None):
        """
        Find a secondary result/attribute from a GrowthModel. Available
        results/attributes can be found in the GrowthModel.outputs() method.

        Parameters
        ----------
        alt_requires : dict, optional
            Alternative dictionary to search for, if you are just interested in
            the oupupts of a GrowthModel child class (for example, cs_outputs
            from ``BioenergeticGrowthModel``). If None, self.requires is used.
            Default None.
        """
        this_req = self.requires
        if alt_requires:
            this_req = alt_requires

        subparams = {}
        for n, loc in this_req.items():
            if isinstance(loc, list):
                for h in loc:
                    if n in h.outputs():
                        subparams[n] = h.outputs()[n]
            elif isinstance(loc, dict):
                for k, h in loc.items():
                    if n in h.outputs():
                        subparams[n] = h.outputs()[n]
            else:
                if n in h.outputs:
                    subparams[n] = h.outputs[n]
        return subparams
