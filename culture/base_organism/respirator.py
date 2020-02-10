"""

This is the respiration submodule, designed to support organism by perfoming
calculations related to its respiration.

@author P M Higgins
@version 0.0.1

"""

import sys
sys.path.append("../..")
import math
from NutMEG import reaction as rxn
from NutMEG.environment import environment
from NutMEG.reactor import reactor

import NutMEG.util.NutMEGparams as nmp
from NutMEG.util.loggersetup import loggersetup as logset
logger = logset.get_logger(__name__, filelevel=nmp.filelevel, printlevel=nmp.printlevel)

class respirator:
    """
    Method class for implementing Jin and Bethke (2007)'s procedure for
    estimating the rate of nutrient and energy uptake. Can be extended for
    specific circumstances e.g. methanogenesis.
    """
    net_pathway = None # reaction describing the net catabolic reaction
    name='pathway'
    ATP_production = None # reaction describing production of ATP
    # locale = reactor() # local chemical environment in which the reaction
      # is taking place. Can be passed as inside or outside the cell for
      # various mechanisms.

    G_A = 0.0 # total free energy of the overall catabolic pathway
      # (per molar overall reaction)
    G_P = 0. # total free energy of each ATP producion
      # (per mol of ATP produced) (will be +ve)
    G_C = 0. # total free energy to be conserved by catbolism
      # (per) molar overall reaction
    n_P = 0.0 # relative total number of ATP formed per pathway
    n_HR = 0.0 # realative total number of +ve ions transferred across membrane
      # per pathway
    n_HP = 3.0 # relative total number of H+ ions translocated per
      # ATP synthesis (usually 3)
    n_ATP = 0.0 # total number of ATP produced per pathways based on
      # the above 3 n's.
    #RTP = 0.035/3600. # rate constant of rds in pathway (usually ATP synthesis)
      # this number is from Moestedt:2015 for methanogenesis, there may be
      # better sources out there
    #k_T = None # rate constant in the environment, no P dependence in yet.
    xi = 1.0 # Stoichiometric coefficient. Avg no of times the rds has taken
      # place. Usually 1.
    F_T = None # scaling factor due to thermodynamic effects.
    max_rate = None # calclated thermodynamically limited rate of reaction
      # in the forwards direction.
    rate = None # the actual rate, corrected for other limiters in
      # the organism.

    def __init__(self, host, net_pathway, n_ATP,
      celldata=[0.0001, 0.004, 0.005, 7.], name='pathway',
      xi=1.0, G_net_pathway=None, pathwaytype=None,
      *args, **kwargs):
        self.host = host
        self.locale = host.locale
        self.name = name
        # unify the pathway with the local environment
        if type(net_pathway) is str:
            if pathwaytype == None:
                pathwaytype = type(rxn.reaction({},{}, self.locale.env))
            self.net_pathway = self.locale.reactionlist[net_pathway][pathwaytype]
        elif type(net_pathway) is rxn.reaction or rxn.redox:
            # add reaction direct to the reactor, if it isn't there already
            self.locale.add_reaction(net_pathway, overwrite=kwargs.pop('overwrite', False))
            #set self.net_pathway now it has been unified
            self.net_pathway = self.locale.reactionlist[net_pathway.equation][type(net_pathway)]
        else:
            raise ValueError('Unable to process your reaction type')
        self.xi=xi
        if G_net_pathway is not None:
            self.G_A = G_net_pathway
        else:
            # get the molar gibbs from the net pathway reaction ourselves
            self.net_pathway.rto_current_env()
            self.net_pathway.update_molar_gibbs_from_quotient(
              updatestdGibbs=False)
            self.G_A = self.net_pathway.molar_gibbs

        self.build_ATP_reaction(celldata) # also gets G_P

        if n_ATP is None:
            self.n_P = kwargs.pop('n_P', 0.0)
            self.n_HP = kwargs.pop('n_HP', 3.0)
            self.n_HR = kwargs.pop('n_HR', 0.0)
            self.n_ATP = (self.n_P + (self.n_HR/self.n_HP))
        else:
            self.n_ATP = n_ATP

        # self.k_RTP = kwargs.pop('k_RTP', 0.035/3600.)
        if self.net_pathway.rate_constant_RTP is None:
            self.net_pathway.rate_constant_RTP = kwargs.pop('k_RTP', 0.0001586)

        if not self.net_pathway.bool_rate_constants():
            raise ValueError('Rate constant unknown, we cannnot' + \
            'calucate respiration rates without them!')

        if self.net_pathway.rate_constant_env == None:
            #if we don't know the rate constant outside RTP,
            # it often changes by 2x every increase by 10 K.
            self.net_pathway.rate_constant_env = ( \
              self.net_pathway.rate_constant_RTP * \
              (2**((self.locale.env.T-298)/10)))

        #self.k_RTP = kwargs.pop('k_RTP', 0.065/60.)
        # self.k_T = self.k_RTP*math.exp(298.0/self.locale.env.T)

        # self.k_T = 34.7/3600#self.k_RTP*math.exp(6000*((1/298)-(1/self.locale.env.T)))
          # ^ get the rate const at this temperature
          # (only accurate for near RTP)

        self.G_C = self.n_ATP*self.G_P


    def build_ATP_reaction(self, celldata):
        """Create a reaction object describing the formation of ATP using
        cell parameters.

        celldata is in the form [activity ADP, activity P, activity ATP, pH]
        """
        ADP = rxn.reagent('+H3ADP1-(aq)', self.locale.env, activity=celldata[0],
          phase='aq')
        P = rxn.reagent('H3PO4(aq)', self.locale.env, activity=celldata[1],
          phase='aq')
        ATP = rxn.reagent('+H4ATP-(aq)', self.locale.env, activity=celldata[2],
          phase='aq')
        #H = reaction.reagent('H+', self.env, activity=(10**-celldata[3]),
        #  phase='aq', molar_ratio=2.)
        H2O = rxn.reagent('H2O(l)', self.locale.env, phase='l',
          conc=55.5, phase_ss=True, activity=1.0)

        self.ATP_production = rxn.reaction({ADP:1, P:1},
          {ATP:1, H2O:1}, self.locale.env)
        self.ATP_production.rto_current_env()

        self.ATP_production.update_molar_gibbs_from_quotient(
          updatestdGibbs=False)
        self.G_P = self.ATP_production.molar_gibbs



    def get_rate(self):
        """Update the rate of reaction in the last state this object
        was left in.
        """
        self.G_A = self.net_pathway.molar_gibbs
        f = (-self.G_A-self.G_C)
        # ^ the net of energy released, thermodynamic driving force.
        if f <= 0:
            self.F_T = 0. # reaction cannot proceed. Technically, it goes backwards
        else:
            self.F_T = 1-math.exp(-f/(self.xi*8.314472*self.locale.env.T))

        conc_multiplier = 1.0
        for r, mr in self.net_pathway.reactants.items():
            conc_multiplier = conc_multiplier*(r.activity**mr)

        self.rate = (self.net_pathway.rate_constant_env * \
          conc_multiplier*self.F_T)
