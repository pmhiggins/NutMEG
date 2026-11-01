from collections.abc import Iterable
import numpy as np

class RateAggregator:
    """
    Superclass for all rate aggregation calculations. All instances must have
    a version of the function ``combine(base_rate, factors)``.
    """

    def combine(self, base_rate, factors):
        """ Calculate and return the aggregated rate."""
        raise NotImplementedError

    @staticmethod
    def unify_iterables(factors):
        """
        Check if at least one of the factors passed is an iterable. If so,
        transform all other factors to have the same dimension as it does
        and return as a numpy array.
        """

        iters = [isinstance(_f, Iterable) for _f in factors]
        if any(iters):
            _n = None
            for _i, _it in enumerate(iters):
                if _it:
                    _n = len(factors[_i])
                    break

            for _i, _f in enumerate(factors):
                if not iters[_i]:
                    factors[_i] = np.array([_f,]*_n)

        return np.array(factors)
