from .rate_aggregator import RateAggregator

class Multiplicative(RateAggregator):
    def combine(self, base_rate, factors):
        r = base_rate
        for f in factors:
            r *= f
        return r
