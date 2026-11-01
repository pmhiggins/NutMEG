from .forcing_factor import ForcingFactor
import numpy as np

class BiologicalPerformance(ForcingFactor):
    r"""
    A forcing factor for implementing a biological performance curve as
    described by Yin et al (1995) and Mendez et al (2021):

    .. math::

        F = \left\[\left\(frac{x-x_{min}}{x_{opt}-x_{min}}\right\)\left\(frac{x_{max}-x}{x_{max}-x_{opt}}\right\)^{\frac{x_{max}-x_{opt}}{x_{opt}-x_{min}}}\right\]^{s}

    where :math:`x` is the controlling parameter (e.g., temperature),
    :math:`x_{opt}` is its optimum value (i.e., when growth is largest), and
    :math:`x_{min}` and :math:`x_{max}` are the minimum and maximum values of x
    for which rate is >0. The parameter :math:`s` controls the shape of the
    resulting bell curve.

    Attributes
    ----------
    x : str
        Reactor attribute that acts as the controlling parameter. options include:
        'T' for temperature, 'P' for pressure, 'pH', or any string identifier
        for any Reagent object.
    x_opt : float
        The optimum value of x
    x_min : float
        The minimum, or 'cutoff' value of x, below which the rate should be 0.
    x_max : float
        The maximum, or 'cutoff' value of x, above which the rate should be 0.
    s : float
        Value of the shape parameter.
    requires : dict or Nonetype
        Dictionary of additional required host properties to run
        this GrowthModel. Keys are property identifiers, and values are the
        object in a NutMEG.core class to look in.
    """

    def __init__(self, x, x_opt, x_min, x_max, s, conctype=None):
        """
        Extends ForcingFactor.__init__()

        Parameters
        ----------
        x : str
            Reactor attribute that acts as the controlling parameter. options include:
            'T' for temperature, 'P' for pressure, 'pH', or any string identifier
            for any Reagent object.
        x_opt : float
            The optimum value of x
        x_min : float
            The minimum, or 'cutoff' value of x, below which the rate should be 0.
        x_max : float
            The maximum, or 'cutoff' value of x, above which the rate should be 0.
        s : float
            Value of the shape parameter.
        conctype : str
            If x is a chemical species, the nature of concentration to use. Can
            be either 'conc', 'molal', or 'activity'. Activity will be used if
            not specified.
        """
        super().__init__()
        self.x_str = x
        self.x_min = x_min
        self.x_max = x_max
        self.x_opt = x_opt
        self.s = s
        self.conctype = conctype
        if self.conctype == None:
            self.conctype = 'molality'



    def get_x_val(self, host, locale):
        if self.x_str =='T':
            return locale.T
        elif self.x_str == 'P':
            return locale.P
        elif self.x_str == 'pH':
            return locale.pH
        elif self.x_str in locale.composition.keys():
            return locale.composition[self.x_str].get_amount(self.conctype)


    def compute(self, host, locale):
        """
        Calculate and return the inhibition forcing factor.

        Parameters
        ----------
        host : base_organism
            Host organism. Can be passed as None for this ForcingFactor if using
            a factor-specific max rate.
        locale : reactor
            Host chemical reactor. If photosynthetically active radiation is
            not set manually to this InhibitionJassbyPlatt, it should be
            up-to-date in locale via `locale.PAR`
        """
        x = self.get_x_val(host, locale)

        if x < self.x_min or x>self.x_max:
            return 0.

        _a = (x-self.x_min)/(self.x_opt - self.x_min)
        _b = (self.x_max-x)/(self.x_max-self.x_opt)
        _b = _b**((self.x_max - self.x_opt)/(self.x_opt - self.x_min))
        F = (_a*_b)**self.s

        return F
