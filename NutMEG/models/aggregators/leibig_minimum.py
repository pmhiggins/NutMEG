from .rate_aggregator import RateAggregator
import numpy as np



class LeibigMinimum(RateAggregator):
    r"""
    RateAggregator that computes a rate and forcing factors according to
    Leibig's law of the minimum, i.e.,:

    .. math::

        r = r_{max} \times \min\left\{F_1, F_2, F_3, ...\right\}

    """

    def combine(self, base_rate, factors):
        """
        Implement Leibig's law of the minimum.

        Parameters
        ----------
        base_rate : float
            Maximum rate
        factors : iterable
            List of outputs from a ``ForcingFactor.compute()`` calculation.
            Each entry may be a scalar or vector (e.g., a 1D array), but if
            multiple vectors are in ``factors``, they must have the same length.

        Returns
        -------
        Aggregated rate as a ``float`` or ``numpy.array`` if element(s) of
        ``factors`` are iterable
        """
        factors_arr = RateAggregator.unify_iterables(factors)
        return base_rate * np.min(factors_arr, axis=0)
