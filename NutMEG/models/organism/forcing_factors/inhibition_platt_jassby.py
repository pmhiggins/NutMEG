from .forcing_factor import ForcingFactor
import numpy as np

class InhibitionPlattJassby(ForcingFactor):
    r"""
    A forcing factor for phototrophy using the model presented by Platt & Jassby
    (1976):

    .. math::

        F = \tanh(\alpha I  /\mu_{max})

    where :math:`I` is the photosynthetically active radiation in
    :math:`\mu mol photons m^{-2} s^{-1}`, `\mu_{max}` is the maximum (growth)
    rate in :math:`s^{-1}` and :math:`\alpha` is a fitting parameter.

    Attributes
    ----------
    alpha : float
        Platt-Jassby fitting coefficient
    mu_max : float or Nonetype
        The value of max rate to be used for this calculations. If None, the
        organism-level maximum growth rate will be used. Default None.
    I : float or Nonetype
        Local flux of photosynthetically active radiation. If None,
        InhibitionJassbyPlatt will attempt to use a `PAR` attribute
        of the local reactor.
    requires : dict or Nonetype
        Dictionary of additional required host properties to run
        this GrowthModel. Keys are property identifiers, and values are the
        object in a NutMEG.core class to look in.
    """

    def __init__(self, alpha, mu_max=None, I=None):
        """
        Extends ForcingFactor.__init__()

        Parameters
        ----------
        alpha : float
            Platt-Jassby fitting coefficient
        mu_max : float, optional
            The value of max rate to be used for this calculations. If None, the
            organism-level maximum growth rate will be used. Default None.
        I : float, optional
            Local flux of photosynthetically active radiation. If None,
            InhibitionJassbyPlatt will attempt to use a `PAR` attribute
            of the local reactor.
        """
        super().__init__()
        self.alpha = alpha
        self.mu_max = mu_max
        self.I = I

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
        _I = self.I
        if self.I == None:
            try:
                _I = locale.PAR
            except:
                raise ValueError('Unable to retrieve PAR. Did you set reactor.PAR?')

        if self.mu_max:
            return np.tanh(self.alpha * _I  /self.mu_max)
        else:
            return np.tanh(self.alpha * _I  /host.growth.max_growth_rate)
