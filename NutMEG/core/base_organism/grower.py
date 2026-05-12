
from NutMEG.models.organism.base_rate_models.constant import Constant
from NutMEG.models.organism.base_rate_models.base_rate_model import BaseRateModel
from NutMEG.models.organism.growth_models.bioenergetic_growth_model import BioenergeticGrowthModel
from NutMEG.models.aggregators.multiplicative import Multiplicative


class Grower:
    """Class for managing organismic growth rates.

    Attributes
    ----------
    base_rate : BaseRateModel
    forcing_factors : list[ForcingFactor]
    aggregator : Aggregator
    growth_model : GrowthModel

    Notes
    -----
    In NutMEG v1, the predecessor to this class also managed nutrient uptake
    from the environment. Grower doesn't do so at the moment, but might in
    the future.

    """


    def __init__(self,
      base_rate = 'default',
      forcing_factors = 'default',
      aggregator = 'default',
      growth_model = 'default',
      forcing_factor_labels = [],
      nutrient_sources={}):

        # self.host = host
        self.base_rate = base_rate
        self.forcing_factors = forcing_factors
        self.aggregator = aggregator
        self.growth_model = growth_model

        if self.base_rate == 'default':
            self.base_rate = BaseRateModel()
        elif self.base_rate is float:
            self.base_rate = Constant(self.base_rate)

        if self.forcing_factors == 'default':
            self.forcing_factors = []
        if self.aggregator == 'default':
            self.aggregator = Multiplicative()
        if self.growth_model == 'default':
            self.growth_model = BioenergeticGrowthModel()


        self.growth_rate = None
        self.max_growth_rate = None
        self.gross_growth_rate = None

        if not forcing_factor_labels:
            self.forcing_factor_labels = range(len(self.forcing_factors))
        else:
            self.forcing_factor_labels = forcing_factor_labels



    def set_max_rate(self, _rate):
        """
        Set the maximum allowable growth rate to _rate, which can be a constant
        value or a BaseRateModel instance.
        """
        if _rate == 'default':
            self.base_rate = BaseRateModel()
        elif type(_rate) is float:
            self.base_rate = Constant(_rate)
        elif isinstance(_rate, BaseRateModel):
            self.base_rate = _rate
        else:
            raise ValueError('Unknown growth rate update passed')

    def compute_rate(self, host, locale):
        """ Use the growth model to estimate the cell-specific growth rate of host."""
        g_out = self.growth_model.compute(host, locale)
        return g_out

    def get_forcing_factors(self, host, locale):
        return {k:f.compute(host, locale) for k,f in zip(self.forcing_factor_labels, self.forcing_factors)}
