from .growth_model import GrowthModel


class ClassicGrowthModel(GrowthModel):
    """
    Superclass for all growth rate calculations.

    Attributes
    ----------
    requires : list
        List of additional required host properties to run
        this GrowthModel.
    """

    def __init__(self, **kwargs):
        super().__init__()


    def compute(self, host, locale):
        """
        Calculate and return the growth rate according to this classic model.
        Does not intrinsically depend on metabolic rate or maintenance, only
        the max growth rate and growth forcing factors.

        Parameters
        ----------
        host : BaseOrganism
            Host organism. must be capable of computing a max_growth_rate.
        locale : Reactor
            Host chemical reactor. Some ForcingFactors will need this to initialise and
            some won't. It is best to assume they will (else they may throw an error)
        """
        host.growth.max_growth_rate = host.growth.base_rate.compute(host, locale)

        values = [f.compute(host, locale) for _,f in host.growth.forcing_factors.items()]
        host.growth.growth_rate =  host.growth.aggregator.combine(host.growth.max_growth_rate, values)
        return host.growth.growth_rate, host.growth.growth_rate
