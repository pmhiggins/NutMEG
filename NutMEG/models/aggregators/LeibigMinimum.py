from .RateAggregator import RateAggregator

class LeibigMinimum(RateAggregator):
    def combine(self, base_rate, factors):
        return base_rate * min(factors)
