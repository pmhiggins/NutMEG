import math
from ForcingFactor import ForcingFactor
from NutMEG import reaction as rxn

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

    def __init__(self, host,
      xi=1.,
      n_ATP=1.,
      G_ATP='default',
      G_C= None,
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
            reaction
        n_ATP : float, optional
            Number of moles of ATP yielded per mole of ``net_pathway``. Default 1.
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


        super().__init__()
        self.xi = xi

        # set the free energy of the ATP synthesis reaction
        if G_ATP == 'default':
            # use deafult: ATP production at RTP with default celldata
            self.G_ATP = 59623.7
        elif type(G_ATP) == type(0.) or type(G_ATP) == type(0):
            # Hard-coded value has been passed.
            self.G_ATP = float(G_ATP)
        if G_ATP is 'build':
            # built an ATP reaction in the host's locale.
            ATP_rxn = self.build_ATP_reaction(celldata)
            self.G_ATP = ATP_rxn.molar_gibbs
        else:
            raise TypeError('Unknown type of G_ATP passed: '+str(type(G_ATP)))

        # set the number of ATP yielded per host.net_pathway
        if n_ATP is None:
            if kwargs.get('n_P',None) and kwargs.get(n_HP, None) and kwargs.get(n_HR, None):
                self.n_ATP = Bioenergetic.get_nATP_from_protons(n_P,n_HP,n_HR)
            elif G_C:
                # conservable gibbs has been passed directly.
                # use this to set a proxy n_ATP
                self.n_ATP = G_C / self.G_ATP
            else:
                raise ValueError('No n_ATP or G_C passed to Bioenergetics ForcingFactor')
        else:
            self.n_ATP = n_ATP

        # set the conservable energy per mole of host.net_pathway
        self.G_C = self.n_ATP * self.G_ATP


    def compute(self, host):
        """
        Overrides ForcingFactor.compute(). Returns the fractional reduction
        in metabolic rate owing to thermodynamic forcing.
        """
        fT = resp.f_T()
        _f = -host.metabolism.G_A - self.G_C
        if _f <=0.:
            return 0.
        else:
            return max(0., 1-math.exp(-(fT)/(self.xi*8.314472*host.locale.T)))


    def outputs(self, host):
        host.respiration.G_C = self.G_C
        return {'host.respiration.G_C':self.G_C}


    @staticmethod
    def get_nATP_from_protons(n_P, n_HR, n_HP):
        return n_P + (n_HR/n_HP)


    def build_ATP_reaction(self, celldata):
        """
        Create a reaction object describing the formation of ATP using
        cell parameters.

        celldata is in the form [activity ADP, activity P, activity ATP, pH]
        """
        if celldata == 'Build':
            celldata = [0.0001, 0.004, 0.005, 7.]

        ADP = rxn.reagent('+H3(ADP)(aq)', host.locale, activity=celldata[0],
          phase='aq')
        P = rxn.reagent('H3PO4(aq)', host.locale, activity=celldata[1],
          phase='aq')
        ATP = rxn.reagent('+H4(ATP)(aq)', host.locale, activity=celldata[2],
          phase='aq')
        #H = reaction.reagent('H+', self.env, activity=(10**-celldata[3]),
        #  phase='aq', molar_ratio=2.)
        H2O = rxn.reagent('H2O(aq)', host.locale, phase='l',
          conc=55.5, phase_ss=True, activity=1.0)

        ATP_production = rxn.reaction({ADP:1, P:1},
          {ATP:1, H2O:1}, self.locale.env)
        ATP_production.rto_current_env()

        ATP_production.update_molar_gibbs_from_quotient(
          updatestdGibbs=False)
        return ATP_production
