from .forcing_factor import ForcingFactor
from collections.abc import Iterable


class LowerLimit(ForcingFactor):
    """
    Step function like forcing factor for microbial kinetic models. Above some
    limit, the factor value is 0., below it is 1.

    Attributes
    ----------
    species : str
        Name of chemical species which is rate-limiting.
    val : float
        Lower limit concentration of that species
    conctype : str
        Concentration identifier of ``reagent`` to use. Can be either 'conc',
        'molal', 'activity'.
    """

    def __init__(self, species, val, conctype='molality'):
        """
        Extends ForcingFactor.__init__()

        Parameters
        ----------
        species : str
            Name of chemical species which is rate-limiting.
        val : float
            Lower limit concentration of that species
        conctype : str
            Concentration identifier of ``reagent`` to use. Can be either 'conc',
            'molal', 'activity'.
        """
        self.species = species
        self.val = val
        self.conctype=conctype


    def compute(self, host, locale):
        """
        Calculate and return the LowerLimit forcing factor.

        Parameters
        ----------
        host : base_organism
            Host organism. Can be passed as None for this ForcingFactor.
        locale : reactor
            Host chemical reactor. Must have a correct concentration of
            ``substrate``.
        """
        S = locale.composition[self.species].get_amount(self.conctype)

        if isinstance(self.val, Iterable):
            _ret = deepcopy(self.val)
            for i in range(len(self.val)):
                _ret[i] = 0. if S < self.val[i] else 1.
            return _ret
        return 0. if S < self.val else 1.
