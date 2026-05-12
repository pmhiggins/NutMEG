
import sys, warnings, math
from NutMEG.core.reactor.reaction import Reaction
from NutMEG.core.reactor.reagent import Reagent

from NutMEG.models.organism.base_rate_models.constant import Constant
from NutMEG.models.organism.base_rate_models.base_rate_model import BaseRateModel
from NutMEG.models.organism.forcing_factors.bioenergetic import Bioenergetic
from NutMEG.models.aggregators.multiplicative import Multiplicative


class Metaboliser:
    """
    Class for handling kinetics of a net metabolic reaction in an organism.

    Attributes
    ----------
    net_pathway : Reaction
        The overall metabolic reaction to use.
    base_rate : BaseRateModel
        Model to be used for calculating the maximum metabolic rate
        in environmental conditions (i.e., before forcing, if present)
    forcing_factors : list[ForcingFactor]
        List of ForcingFactor models to use to estimate the actual metabolic
        rate.
    aggregator : RateAggregator
        method to use to select the actual rate from the forcing factors
        (e.g., Multiplicative, LeibigMinimum, etc.)
    max_rate : float
        the maximum possible reaction rate for the overal metabolic reaction.
    rate : float
        the actual reaction rate, corrected for other limiters in
        the organism.
    """

    def __init__(self, net_pathway,
      base_rate = 'default',
      forcing_factors = 'default',
      aggregator = 'default',
      forcing_factor_labels = None,
      overwrite_net_pathway=False):
        """
        Parameters
        ----------
        host : BaseOrganism
            host organism.
        locale : Reactor
            The chemical reactor the organism exists inside.
        net_pathway : Reaction or str
            The overall metabolism to use. In any case it is best to pass a
            reaction like object (reaction or redox) which will then be unified
            with ``locale``. If a string is passed, look in ``host.locale``
            for the reaction and set that as the pathway.
        base_rate : BaseRateModel
            Model to be used for calculating the maximum metabolic rate
            in environmental conditions (i.e., before forcing, if present)
        forcing_factors : list[ForcingFactor]
            List of ForcingFactor models to use to estimate the actual metabolic
            rate.
        aggregator : RateAggregator
            method to use to select the actual rate from the forcing factors
            (e.g., Multiplicative, LeibigMinimum, etc.)
        forcing_factor_labels : list, optional
            List of identifiers for forcing_factors.
        overwrite_net_pathway : bool, optional
            Pass if net_pathway has been created but not yet unified with
            the host's locale.
        """

        self.net_pathway = net_pathway

        #### set up kinetics
        self.base_rate = base_rate
        self.forcing_factors = forcing_factors
        self.aggregator = aggregator

        if self.base_rate == 'default':
            self.base_rate = BaseRateModel()
        elif self.base_rate is float:
            self.base_rate = Constant(self.base_rate)
        if self.forcing_factors == 'default':
            self.forcing_factors = [Bioenergetic()]
        if self.aggregator == 'default':
            self.aggregator = Multiplicative()


        if not forcing_factor_labels:
            self.forcing_factor_labels = range(len(self.forcing_factors))
        else:
            self.forcing_factor_labels = forcing_factor_labels

        self.rate = None
        self.max_rate = None


    def set_max_rate(self, _rate):
        if _rate == 'default':
            self.base_rate = BaseRateModel()
        elif type(_rate) is float:
            self.base_rate = Constant(_rate)
        elif isinstance(_rate, BaseRateModel):
            self.base_rate = _rate
        else:
            raise ValueError('Unknown metabolic rate update passed')

    def compute_rate(self, host, locale):
        """
        Calulate the current maximum metabolic rate. This needs to be
        computed every time because some max rate methods depend on the current
        environmental conditions.
        """
        self.max_rate = self.base_rate.compute(host, locale)

        values = [f.compute(host, locale) for f in self.forcing_factors]
        self.rate =  self.aggregator.combine(self.max_rate, values)
        return self.rate

    def get_forcing_factors(self, host, locale):
        return {k:f.compute(host, locale) for k,f in zip(self.forcing_factor_labels, self.forcing_factors)}
