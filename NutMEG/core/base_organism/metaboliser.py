
import sys, warnings
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
    Class for implementing respiration in an organism. Uses Jin and Bethke
    (2007)'s procedure for estimating the rate of nutrient and energy uptake.
    Can be extended for specific circumstances e.g. methanogenesis.

    Attributes
    ----------
    host : ``base_organism`` like
        host organism. Ensure that the host organism's locale object is the
        reactor you want.
    net_pathway : ``reaction`` like or str
        The overall metabolism to use. In any case it is best to pass a
        reaction like object (reaction or redox) which will then be unified
        with ``host.locale``. If a string is passed, look in ``host.locale``
        for the reaction and set that as the pathway.

    name : str, optional
        Name of the pathway. Default is 'pathway', used for selecting the
        reaction from ``host.locale`` if ``net_pathway`` is passed as a str.

    ATP_production : ```reaction`` like
        reaction from the porduction of ATP.
    G_A : float
        total free energy of the overall catabolic pathway (per molar overall
        reaction)
    F_T : float
        Scaling factor due to thermodynamic effects.
    rate : float
        the actual reaction rate, corrected for other limiters in
        the organism.


    """
    # net_pathway = None # reaction describing the net catabolic reaction
    # name='pathway'
    # ATP_production = None # reaction describing production of ATP
    # locale = reactor() # local chemical environment in which the reaction
      # is taking place. Can be passed as inside or outside the cell for
      # various mechanisms.

    # G_A = 0.0 # total free energy of the overall catabolic pathway
    #   # (per molar overall reaction)
    # G_P = 0. # total free energy of each ATP producion
    #   # (per mol of ATP produced) (will be +ve)
    # G_C = 0. # total free energy to be conserved by catbolism
    #   # (per) molar overall reaction
    # n_P = 0.0 # relative total number of ATP formed per pathway
    # n_HR = 0.0 # realative total number of +ve ions transferred across membrane
    #   # per pathway
    # n_HP = 3.0 # relative total number of H+ ions translocated per
    #   # ATP synthesis (usually 3)
    # n_ATP = 0.0 # total number of ATP produced per pathways based on
      # the above 3 n's.
    #RTP = 0.035/3600. # rate constant of rds in pathway (usually ATP synthesis)
      # this number is from Moestedt:2015 for methanogenesis, there may be
      # better sources out there
    #k_T = None # rate constant in the environment, no P dependence in yet.
    # xi = 1.0 # Stoichiometric coefficient. Avg no of times the rds has taken
      # place. Usually 1.
    # F_T = None # scaling factor due to thermodynamic effects.
    # max_rate = None # calclated thermodynamically limited rate of reaction
      # in the forwards direction.
    # rate = None # the actual rate, corrected for other limiters in
      # the organism.

    def __init__(self, host, net_pathway,
      n_ATP=1., max_metabolic_rate=None,
      rate_func='first order', rate_func_args={},
      kinetic_forcing_parameters=None, kinetic_F_attrs=None,
      rate_constant_env=None, rate_constant_RTP=None,
      celldata=[0.0001, 0.004, 0.005, 7.], G_ATP=None,
      overwrite_net_pathway=False, G_net_pathway=None, pathwaytype=None,
      n_T=None, n_HP=None, n_HR=None, G_C=None):
        """
        Parameters
        ----------
        celldata : list
            Concentration in the form [activity ADP, activity P, activity ATP, pH]
        G_net_pathway : NoneType or float
            If you want to pass the gibbs gree energy of the metabolic reaction
            use this, if not leave as None. Default ``None``.
        """
        self.host = host
        self.locale = host.locale


        #### unify the pathway with the local environment
        if overwrite_net_pathway:
            if type(net_pathway) is str:
                if pathwaytype == None:
                    pathwaytype = type(rxn.reaction({},{}, self.locale.env))
                self.net_pathway = self.locale.reactionlist[net_pathway][pathwaytype]
            elif type(net_pathway) is rxn.reaction or rxn.redox:
                # add reaction direct to the reactor, if it isn't there already
                self.locale.add_reaction(net_pathway, overwrite=overwrite_net_pathway)
                #set self.net_pathway now it has been unified
                self.net_pathway = self.locale.reactionlist[net_pathway.equation][type(net_pathway)]
            else:
                raise ValueError('Reactor of ',self.host.name, 'is unable to process net_pathway type')
        else:
            self.net_pathway = self.locale.reactionlist[net_pathway.equation][type(net_pathway)]

        #### set up overall free energy of metabolic reaction

        if G_net_pathway:
            self.G_A = G_net_pathway
            self.net_pathway.molar_gibbs = self.G_A
        else:
            # get the molar gibbs from the net pathway reaction ourselves
            self.net_pathway.rto_current_env()
            self.net_pathway.update_molar_gibbs_from_quotient(
              updatestdGibbs=False)
            self.G_A = self.net_pathway.molar_gibbs



        #### set up metabolic rates

        # if a rate constant is passed, the net_pathway's rate constant will be
        # overwritten. If not, it will be left untouched.
        if rate_constant_env:
            self.net_pathway.rate_constant_env=rate_constant_env
        if rate_constant_RTP:
            self.net_pathway.rate_constant_RTP=rate_constant_RTP

        if not self.net_pathway.bool_rate_constants():
            warnings.warn('Respirator initiated without any rate constants')
             # unknown, we cannnot' + \
            # 'calucate respiration rates without them!')
        else:
            if self.net_pathway.rate_constant_RTP and not self.net_pathway.rate_constant_env:
                # if we don't know the rate constant outside RTP,
                # it often changes by 2x every increase by 10 K.
                self.net_pathway.rate_constant_env = ( \
                  self.net_pathway.rate_constant_RTP * \
                  (2**((self.locale.env.T-298)/10)))

        #### set up rate function.
        self.max_metabolic_rate = max_metabolic_rate
        if rate_func == 'first order':
            self.rate_func_ID = 'first order'
            self.rate_func = self._rf_first_order
        elif rate_func == 'Arrhenius':
            self.rate_func_ID = 'Arrhenius'
            self.rate_func = self.net_pathway.calculate_rate
        elif rate_func == 'zeroth order':
            self.rate_func_ID = 'zeroth order'
            self.rate_func = lambda: self.net_pathway.rate_constant_env
            self.max_metabolic_rate = self.net_pathway.rate_constant_env
        else:
            self.rate_func = rate_func
        self.rate_func_args = rate_func_args


        #### setup forcing parameters for respiration
        # Load default forcing functions
        self.forcing_parameters = {}
        self.F_attrs = {}
        if not kinetic_F_attrs:
            self.F_attrs = {'xi':1.0}
        else:
            self.F_attrs = kinetic_F_attrs
            self.F_attrs['xi'] = kinetic_F_attrs.get('xi', 1.0)
        self._set_default_forcing()

        if kinetic_forcing_parameters:
            for name, (func, arg_keys) in kinetic_forcing_parameters.items():
                self.set_forcing_parameter(name, func, arg_keys)

        self.get_rate()








    def set_forcing_parameter(self, name, func, arg_keys):
        """
        Set or override a forcing function and define its required arguments.
        Note that your function MUST take a respiration object as its first argument.
        Names of remaining arguments should be listed in arg_keys. Be sure to update
        respirator.F_attrs with the arguments you wish to use, or you will recieve a cryptic error.
        """
        if not callable(func):
            raise ValueError(f"Forcing function for {name} must be callable.")
        if not isinstance(arg_keys, list):
            raise TypeError(f"Argument keys for {name} must be a list of parameter names.")

        self.forcing_parameters[name] = (func, arg_keys)




    def get_forcing_fraction(self, F_ID):
        """
        For forcing parameter with identifier F_ID, return the instantaneous
        rate forcing (usually between 0 and 1).
        """
        try:
            func, arg_keys = self.forcing_parameters[F_ID]
            args = [self.F_attrs[key] for key in arg_keys if key in self.F_attrs]
            return func(self, *args)  # Apply the forcing function
        except:
            warnings.warn("Forcing parameter: '"+F_ID+" not present for "+self.host.name)
            return None


    def get_rate(self):
        """
        Update the rate of reaction in the last state this object
        was left in. Applies kinetic forcing parameters based on environement,
        if they are included in the respirator'a forcing_parameters attribute.
        """
        self.G_A = self.net_pathway.molar_gibbs

        rate_modifier = 1.
        # loop through all defined forcing functions
        # this is implementing the multiplicative monod model.
        for name, (func, arg_keys) in self.forcing_parameters.items():
            # Extract necessary arguments from the environment
            args = [self.F_attrs[key] for key in arg_keys if key in self.F_attrs]
            rate_modifier *= func(self, *args)  # Apply the forcing function
            # print(arg_keys, func(self, *args), rate_modifier)

        self.rate = rate_modifier * self.rate_func(*self.rate_func_args)

        # if there is a hard-coded max rate, check that we have not exceeded it.
        if self.max_metabolic_rate:
            if self.rate > self.max_metabolic_rate:
                self.rate = self.max_metabolic_rate



    def metabolic_energy_density(self):
        """
        Return an approximation of the energy density [J/kg H2O] for the
        net_pathway.

        Calculates the smallest energy yield from 'using up' the metabolic
        reagents. In reality, the free energy would change as the concentration
        decreases, so this is only a measure of the energy density available
        for this metabolism at this moment in time.
        """
        ED = []
        for r, mr in self.net_pathway.reactants.items():
            if r.name != 'H2O(aq)' and r.name != 'H+' and r.name != 'OH-':
                ED.append(r.conc*-self.G_A/mr)
        return min(ED)


    def _rf_first_order(self):
        """
        Generic first order rate law to be added to rate_laws.
        """
        conc_multiplier = 1.0
        for r, mr in self.net_pathway.reactants.items():
            conc_multiplier = conc_multiplier*(r.activity**mr)

        return self.net_pathway.rate_constant_env * conc_multiplier
