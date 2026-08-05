from .forcing_factor import ForcingFactor
from ...aggregators.rate_aggregator import RateAggregator
import numpy as np

class ComboAdd(ForcingFactor):
    r"""
    Compute an additive forcing factor by adding together multiple ForcingFactor
    values. An example use-case is when there are two
    sources of a key nutrient and both can be used used, e.g.:

    .. math::

        F = \min\left\{(F_1 + F_2 + ...), m\right\}

    where m is a hard maximum, default 1.0.

    Parameters
    ----------
    FFs : iterable
        List of ForcingFactor objects to be summed together.
    hard_min : float
        Fixed minimum value this ForcingFactor can return. Default 0.
    hard_max : float
        Fixed minimum value this ForcingFactor can return. Default 1.
    """

    def __init__(self, FFs, hard_min=0., hard_max=1.):
        self.FFs = FFs
        self.hard_min = hard_min
        self.hard_max = hard_max

    def compute(self, host, locale):
        """Calculate and return the value of this ForcingFactor"""
        factors = [_FF.compute(host, locale) for _FF in self.FFs]
        factors = RateAggregator.unify_iterables(factors) # if one factor is iterable, make all have the same dimension
        res = np.sum(factors, axis=0)
        res = np.clip(res, a_min=self.hard_min, a_max=self.hard_max)
        return res
