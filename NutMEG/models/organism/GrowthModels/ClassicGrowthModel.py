from .GrowthModel import GrowthModel


class ClassicGrowthModel(GrowthModel):
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
        super().__init__(host,locale)


    def compute(self, host, locale):
        """ Calculate and return the growth rate according to this model."""
        host.growth.max_growth_rate = host.growth.base_rate.compute(host, locale)

        values = [f.compute(host, locale) for f in host.growth.forcing_factors]
        host.growth.growth_rate =  host.growth.aggregator.combine(host.growth.max_growth_rate, values)
        return host.growth.growth_rate, host.growth.growth_rate
