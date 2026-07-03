from .rate_aggregator import RateAggregator

class Multiplicative(RateAggregator):
    r"""
    RateAggregator that computes a rate and forcing factors according to
    a multiplicative law on all forcing factors i.e.,:
    .. math::
        r = r_{max} \times \prod_{i} F_i

    """

    def combine(self, base_rate, factors):
        """
        Implement the multiplicative rate law.

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
        r = base_rate
        for f in factors:
            r *= f
        return r
