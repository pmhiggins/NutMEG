import math
import warnings
from .ForcingFactor import ForcingFactor
from NutMEG.core.reactor.reaction import reaction as rxn

class Bioenergetic(ForcingFactor):
    """
    Calucates bioenergetic (thermodynamic) forcing factor for microbial
    kinetic models using Jin & Bethke (2007)'s procedure.

    Attributes
    ----------
    xi : float, optional
        Stoichiometric coefficient. The averge no of times the rate-
        determining-step has taken place. Default 1.
    G_ATP : float or str, optional
        Total free energy of each ATP producion (per mol of ATP produced).
        Default 'default'; which sets G_ATP to 59623.7
    G_C : float, optional
        total free energy to be conserved by catbolism per molar overall
        reaction
    n_ATP : float, optional
        Number of moles of ATP yielded per mole of ``net_pathway``. Default 1.
    """

    def __init__(self, host, locale,
      xi=1.,
      n_ATP=None,
      G_ATP='default',
      G_C='default',
      celldata=[0.0001, 0.004, 0.005, 7.], **kwargs):
        """
        Note: At least two of n_ATP, G_ATP, and G_C must be passed for
        successful initialisation.

        Parameters
        ----------
        host : ``base_organism`` like
            host organism. Only required for this ForcingFactor if requesting
            to build the ATP reaction, otherwise can be passed as None.
        xi : float, optional
            Stoichiometric coefficient. The averge no of times the rate-
            determining-step has taken place. Default 1.
        G_ATP : float or str, optional
            Total free energy of each ATP producion (per mol of ATP produced).
            Default 'default'; which sets G_ATP to 59623.7
        G_C : float, optional
            total free energy to be conserved by catbolism per molar overall
            reaction. If passed as a float or int, this will override n_ATP.
            Default behaviour is to take a value of 95% the molar free energy
            of the metabolic yield when the forcing factor is first computed.
        n_ATP : float, optional
            Number of moles of ATP yielded per mole of ``net_pathway``. If both
            n_ATP and G_C are passed as a float, G_C will be prioritised. Default None
        celldata : list, optional
            Concentration in the form [activity ADP, activity P, activity ATP,
            pH]. Only required if requesting to build the ATP reaction (and
            even then, the default should be sufficient in most applications).
        n_P : float, kwarg
            relative total number of ATP formed per pathway Default 0.0
        n_HR : float, kwarg
            realative total number of +ve ions transferred across membrane
            per pathway. Default 0.0
        n_HP : float, kwarg
            relative total number of H+ ions translocated per
            ATP synthesis Default 3.0.
        """


        super().__init__(host, locale)
        self.xi = xi

        # set the free energy of the ATP synthesis reaction
        if G_ATP == 'default':
            # use deafult: ATP production at RTP with default celldata
            self.G_ATP = 59623.7
        elif type(G_ATP) == type(0.) or type(G_ATP) == type(0):
            # Hard-coded value has been passed.
            self.G_ATP = float(G_ATP)
        elif G_ATP == 'build':
            # built an ATP reaction in the host's locale.
            ATP_rxn = self.build_ATP_reaction(locale, celldata)
            self.G_ATP = ATP_rxn.molar_gibbs
        else:
            raise TypeError('Unknown type of G_ATP passed: '+str(type(G_ATP)))


        self.default_GC_init = False

        if n_ATP is None:
            if type(G_C) is float or type(G_C) is int:
                self.G_C = float(G_C)
                self.n_ATP = self.G_C / self.G_ATP
            elif G_C == 'default':
                # somehow set G_C and n_ATP based on net_metabolism. Add a catch
                # for compute so it does so then, after resp is initialised?
                self.default_GC_init = True
            elif kwargs.get('n_P',None) and kwargs.get(n_HP, None) and kwargs.get(n_HR, None):
                self.n_ATP = Bioenergetic.get_nATP_from_protons(n_P,n_HP,n_HR)
                self.G_C = self.G_ATP * self.n_ATP
            else:
                raise ValueError('No n_ATP or G_C passed to Bioenergetics ForcingFactor')

        else:
            if type(G_C) is float or type(G_C) is int:
                warnings.warn("Both G_C and n_ATP have been passed to Bioenegetic. Using G_C (n_ATP may change).")
                self.G_C = float(G_C)
                self.n_ATP = self.G_C / self.G_ATP
            else:
                self.n_ATP = float(n_ATP)
                self.G_C = self.G_ATP * self.n_ATP




    def compute(self, host, locale):
        """
        Overrides ForcingFactor.compute(). Returns the fractional reduction
        in metabolic rate owing to thermodynamic forcing.
        """
        if self.default_GC_init:
            # it was not possible to initialise default G_C in __init__, try now
            self.G_C = -0.95*host.metabolism.net_pathway.molar_gibbs
            self.n_ATP = self.G_C / self.G_ATP
            self.default_GC_init == False

        _f = (-host.metabolism.net_pathway.molar_gibbs) - self.G_C
        if _f <=0.:
            return 0.
        else:
            return max(0., 1-math.exp(-(_f)/(self.xi*8.314472*locale.env.T)))


    def outputs(self):
        return {'G_C':self.G_C}


    @staticmethod
    def get_nATP_from_protons(n_P, n_HR, n_HP):
        return n_P + (n_HR/n_HP)


    def build_ATP_reaction(self, locale, celldata):
        """
        Create a reaction object describing the formation of ATP using
        cell parameters.

        celldata is in the form [activity ADP, activity P, activity ATP, pH]
        """
        if celldata == 'Build':
            celldata = [0.0001, 0.004, 0.005, 7.]

        ADP = rxn.reagent('+H3(ADP)(aq)', locale, activity=celldata[0],
          phase='aq')
        P = rxn.reagent('H3PO4(aq)', locale, activity=celldata[1],
          phase='aq')
        ATP = rxn.reagent('+H4(ATP)(aq)', locale, activity=celldata[2],
          phase='aq')
        #H = reaction.reagent('H+', self.env, activity=(10**-celldata[3]),
        #  phase='aq', molar_ratio=2.)
        H2O = rxn.reagent('H2O(aq)', locale, phase='l',
          conc=55.5, phase_ss=True, activity=1.0)

        ATP_production = rxn.reaction({ADP:1, P:1},
          {ATP:1, H2O:1}, locale.env)
        ATP_production.rto_current_env()

        ATP_production.update_molar_gibbs_from_quotient(
          updatestdGibbs=False)
        return ATP_production
