"""
This is the solutions submodule of the reaction module. It emphasises
the energy availability of reactions which take place in solutions,
with a few corrections for temperature differences and a degree of
non-ideality. There is a neutralsol class, an electrolyte class and a
redox class. If it gets too messy then I'll switch them up a bit.

All thermodynamic calculations in solution are to be done in molal
form. For now, we are only considering water as the solvent.

@author P M Higgins
@version 0.2
"""

import sys
import os

sys.path.append(os.path.dirname(__file__) + '/..')
sys.path.append(os.path.dirname(__file__) + '/.../environment')
import reaction
import reagent
import numpy as np
import math
from itertools import chain
from warnings import warn
from NutMEG.environment import environment

class neutralsol(reaction.reaction):
    """
    Class for neutral solutes in neutral solutions, and so is rather
    general. The total activity of dissolved solutes can be
    calculated from the osmotic coefficient, which requires the
    activity of the solvent to be known.
    """


    """
        INITIALISATION
    """

    phi = 0. # osmotic coefficient of the solution
    total_molality = 0.
    solvent_activity =1.



    #initialise as a reaction, but only wih the relevant attributes.
    # add other things in as we go along.
    def __init__(self, reactants, products, env, solvent_activity=1.,
      *args, **kwargs):
        reaction.reaction.__init__(self, reactants, products, env,
          *args, **kwargs)

        self.solvent_activity = solvent_activity
        for r in chain(list(self.reactants), list(self.products)):
            if r.phase == 'aq':
                self.total_molality += r.molal
        #calculate the osmotic coefficient
        self.phi = ((55.506*math.log(self.solvent_activity))
          / (-self.total_molality))

    # phi is a function of m so to do this we need to integrate.
      # The solution is polynomial but dependent on our solution.
      # Hence this is not really generally calculable unless we have
      # this polynomial.
      # I'll leave this for now, but if we actually end up really
      # needing it I'll come back.

    def solute_activity(self):
        x = self.phi














class electrolyte(reaction.reaction):
    """
    Class for charged solutes in pH neutral solutions.

    Used to calculate the mean activity/coefficients of
    fully dissociated salts in solution.
    """

    H2O = True # Is the solvent H2O?
    kappa_cons = 2.529122732e20 #constants in the kappa squared expression
    kappa = None # the thickness of the ionic atmosphere
    mean_gamma = None


    def __init__(self, reactants, products, env, H2O=True, *args, **kwargs):
        reaction.reaction.__init__(self, reactants, products, env,
          *args, **kwargs)
        self.H2O = H2O


    def get_permittivity(self):
        """Return the relative pertmittivity of water at temperature
        T, using the expression fitted by Catenacio et al 2003.
        Definitely valid between 250-370 K --- not sure about outside
        that, or at diferent pressures.
        """
        if self.H2O == True:
            return((-0.8292e-6 *(self.env.T**3)) + (0.1417e-2 * (self.env.T**2)) - (0.9297*self.env.T) + (233.76 + 5321./self.env.T))
        else:
            raise ValueError("This is not aqueous and we don't have "
              "the permittivity relationship of your solvent!")

    def get_ionic_strength_c(self):
        """Return the ionic strength of the solution in mol/L.

        Uses the concentration of the reagents.
        """
        I_c = 0.
        for r in chain(list(self.reactants), list(self.products)):
            if r.phase == "aq":
                I_c += ((r.charge**2)*r.conc)
        return(I_c * 0.5)

    def get_ionic_strength_m(self):
        """Return the ionic strength of the solution in mol/kg.

        Uses the molality of the reagents.
        """
        I_m = 0.
        for r in chain(list(self.reactants), list(self.products)):
            if r.phase == "aq":
                I_m += ((r.charge**2)*r.molal)
        return(I_m * 0.5)



    def update_kappa(self):
        """Calculate the kappa of the electrolyte for use in
        Debeye-Huckel theory and MSA theory.

        If we have the volumes, calculations are done using molality
        and density. Otherwise, it uses concentrations.
        """
        epsilon_r = self.get_permittivity()
        if self.env.V == None or self.env.V_RTP == None:
            # we don't have the volume of our container and so can't calculate
              # kappa using our molal properties. Use concentration instead.
            I_c = self.get_ionic_strength_c()
            self.kappa = math.sqrt(self.kappa_cons * 1000. * I_c
              /(epsilon_r * self.env.T)) # where the 1000 is L in a m^-3
        else:
            # calculate the density of water given the container volume V
            # density = density(RTP) * V(RTP) / V(T,P),
              # density(RTP) ~ 1000 for lowish concentrations.
            dens = 1000. * self.env.V_RTP / self.env.V
            I_m = self.get_ionic_strength_m()
            self.kappa = math.sqrt(self.kappa_cons * dens * I_m
              /(epsilon_r * self.env.T)) #where the 1000 is L in a m^-3



    def update_gammas(self):
        """Update the activity coefficent of the aqueous species in
        the reaction using Debeye-Huckel theory.

        Requires reactants to have: conc, phase, charge
        """
        for r in chain(list(self.reactants), list(self.products)):
            if r.conc > 0.05:
                warn('Your concentrations are a little high to be using '
                  'Debeye-Huckel theory... Take your results with a pinch '
                  'of salt.')
        self.update_kappa()
        epsilon_r = self.get_permittivity()
        factor = -self.kappa /(epsilon_r * self.env.T * 3.072356e-33)
        # the number above is 8\pi \epsilon_0 k_B
        for r in chain(list(self.reactants), list(self.products)):
            if r.phase == 'aq' and r.charge != 0:
                r.set_activitycoefficient( math.exp(factor
                  * math.pow((1.602214e-19 * r.charge), 2)))



    def update_mean_gamma(self):
        """Update the mean activity coefficient for the aqueous ions,
        using the mean-spherical approximation.

        Requires reactants to have: conc, radius, charge
        Ideally also pass the volume(s) of the solution to get molal
        activity coefficient.
        """
        for r in chain(list(self.reactants), list(self.products)):
            if r.conc > 1.:
                warn('Your concentrations are a little high to be using '
                  'the MSA... Take your results with more than a pinch '
                  'of salt!')
        sumconc = 0. # sum of concentrations of charged species
        sumrad_p = [] # radii of positively charged species.
        sumrad_n = [] # radii of negatively charged species.
        for r in chain(self.reactants, self.products):
            if r.phase == 'aq' and r.charge != 0:
                if r.radius is None:
                    raise ValueError('You have not provided us with the ionic '
                      'radius of ' + str(r.name))
                sumconc += r.conc
                if r.charge>0:
                    sumrad_p.append(r.radius)
                if r.charge<0:
                    sumrad_n.append(r.radius)
        if self.env.V == None or self.env.V_RTP == None:
            density = sumconc * 6.022e23 * 1000.
            # number density of charged species in m^(-3)
            # assumes the solution has density of 1 kg /L - so RTP only
            # potentially for future use: H2O has 55.506 mol/kg
        else:
            density = sumconc * 6.022e23 * 1000. * self.env.V_RTP / self.env.V
            # because density scales as Volume
        sumrad = np.mean(sumrad_p) + np.mean(sumrad_n)
        # sum of average ionic radii
        #sumrad = 320e-12
        self.update_kappa()
        x = self.kappa * sumrad # the x parameter in MT's derivation
        lng_el = (((x*math.pow((1.+(2.*x)), 0.5)) - x - (x**2))
          / (4.*math.pi*density*(sumrad**3))) # electric contribution
        y = math.pi*density*(sumrad**3)/6.
        lng_HS = ((4.*y) - ((y**2)*9./4.) +((y**3)*3./8.)) / ((1-(y/2.))**3)
        # ^ hard sphere contribution

        self.mean_gamma = math.exp(lng_el + lng_HS)

    def get_a_salt(self, getgamma=True):
        """Return the value of a_{salt}=a_{\pm}^{\nu} for this fully
        dissolved salt.

        a_{\pm} = mean_gamma * mean_molal
        Set getgamma to False if you want to provide your owm mean
        activity coefficient.
        """
        if getgamma == True:
            self.update_mean_gamma()
        m_prod = 1.
        solute_counter= 0
        for r, mr in chain(self.reactants.items(), self.products.items()):
            if r.phase == 'aq' and r.charge != 0:
                m_prod = m_prod * (r.molal**mr)
                solute_counter += mr
        return ((self.mean_gamma**solute_counter) * m_prod)












class redox_half(electrolyte):
    """Class for redox half equations, an extension of electrolyte so it
    can include non-dissociation reactions.

    e.g. of the form Fe3+ + e- => Fe2+
    """

    stdE_RTP = 0.
    stdE = 0.
    dE_by_dT = None
    d2E_by_dT2 = None
    Fcons= 96485.3329
    n = 0. # number of electrons transferred as written
    H2 = None

    def __init__(self, reactants, products, env, n, stdE_RTP, dE_by_dT=None,
      d2E_by_dT2=None, DeltaC_P=0., *args, **kwargs):
        electrolyte.__init__(self, reactants, products, env,
          *args, **kwargs)
        self.stdE_RTP = stdE_RTP
        self.n = n
        self.dE_by_dT = dE_by_dT
        self.d2E_by_dT2 = d2E_by_dT2
        self.DeltaC_P = DeltaC_P

        self.H2 = reagent.reagent('H2(g)', env)

    def update_E(self, estimate=True, newenv=False):
        """Update the electrode potential.

        If we have dE/dT, great.
        If we have d2E/dT2 even better!
        If not, we can estimate them from entropy and isbaric heat
        capacities. Only do so when prompted though...

        ref. salvi and deBethune 1961, Bratsch 1989
        May be worth a recode as using nonstandard temperatures might
        improve accuracy.
        """
        if self.dE_by_dT is not None and estimate==False:
            if self.d2E_by_dT2 is None:
                # calculate the second order derivative from heat capacity
                # zeroes if nothing is passed for heat capacity
                self.d2E_by_dT2 = (self.DeltaC_P
                  /(self.env.T_RTP*self.n*self.Fcons))
            # update
            self.stdE = (self.stdE_RTP
              + ((self.env.T-self.env.T_RTP)*self.dE_by_dT)
              + (0.5*(((self.env.T-self.env.T_RTP)**2)*self.d2E_by_dT2)))
        elif estimate==True:
            if self.thermo and newenv:
                # use the data if we have it and the env has changed.
                # let's update the entropy and Cp of our reagents in the
                  # current environment.
                self.H2.import_params_db()
                for rp in chain(list(self.reactants), list(self.products)):
                    if rp.name != 'e-':
                        rp.import_params_db()
            # we need the enthalpy change of reaction and the heat
              # capacity change of reaction at 298.15 K
            self.update_std_molar_entropy_of_reaction()
            DeltaS = (self.std_molar_entropy
              - (self.n*self.H2.std_formation_entropy_env/2.))
            # ^ take off entropy of H2 (std electrode) (from reaktoro)
            self.dE_by_dT = DeltaS/(self.n*self.Fcons)

            #now the second order...
            self.DeltaC_P = 0.
            for p, mr in self.products.items():
                if p.name != 'e-': #ignore the electron component
                    self.DeltaC_P += (p.Cp_env * mr)
            for r, mr in self.reactants.items():
                if r.name != 'e-':
                    self.DeltaC_P -= (r.Cp_env * mr)
            self.DeltaC_P -= (self.n*self.H2.Cp_env/2.)
            # ^ take off Cp of H2 (std eletrode) (from reaktoro)
            self.d2E_by_dT2 = self.DeltaC_P/(self.env.T*self.n*self.Fcons)

            # finally time to update
            self.stdE = (self.stdE_RTP
              + ((self.env.T-self.env.T_RTP)*self.dE_by_dT)
              + (0.5*(((self.env.T-self.env.T_RTP)**2)*self.d2E_by_dT2)))


        else:
            warn("We don't have a means to calculate the standard electrode "
              "potential at temperatures other than 298K. Using your standard "
              "value --- proceed with caution!")
            self.stdE = self.stdE_RTP




    def update_std_molar_gibbs(self):
        """Update the standard molar gibbs free energy of this half reaction.
        """
        if self.stdE == None:
            raise ValueError("Please first update the electrode potential "
              "with any tabulated data you have before trying to calculate "
              "the free energy!")
        self.std_molar_gibbs = -self.n * self.Fcons * self.stdE








class redox(reaction.reaction):
    """
    Class for redox reactions in solution

    Will make extensive use of electrolye class, so in an ideal world
    our reactants have conc, molal, radius, and the V and V_RTP
    parameters arise from somewhere.

    For now, we are considering a changing temperature and constant
    pressure at 1 bar, for reading off tables, though with reaktoro our
    entropies and heat capacities are updated with pressure. Whether
    this is the whole story is unlikely.

    In the future if we decide to consider buffers or dissociating
    acids etc, these methods can be extended. For now they are limited
    to dissociated salts.

    NOTE: This class only deals in dissociations/half reactions, if you
    want to work out the full redox potentials,use two redox objects
    and find the difference like you would on pen and paper.
    """

    forward = None # reaction.electrolyte object describing the
      # forward reaction solid -> electrolytes
    reverse = None # reaction.electrolyte object describing the
      # reverse reaction solid -> electrolytes
    H2O=True # our solvent
    stdE_RTP = None # standard electrode potential of our couple, in V.
    stdE = None # standard electrode potential at temperature T, in V
    E = None # electrode potential at some non-RTP environment.
    n = 0.0 # number of electrons transferred / reaction
    Fcons= 96485.3329 # Faraday constant C/mol
    env = environment()


    def __init__(self, forward, reverse, cell, n, env=environment(),
      H2O=True, *args, **kwargs):
        # recall that in redox calculations, the forward reaction is the one
          # where it goes aq -> s, i.e. a reaction.electrolyte going
          # backwards! Vice versa for reverse.
        self.forward = forward # forward reaction as a redox half
        self.reverse = reverse # reverse reaction as a redox half
        self.cell = cell    # the full cell reaction, probably as a
          # reaction.reaction object.
          # Useful for clarity, and housing the environment parameters.
          # Still important that the molalities etc are correct though,
          # this is the one we'll be updating!

        self.equation = self.cell.get_equation()
        # assert that each of these has the same surroundings
        self.forward.env = env
        self.reverse.env = env
        self.cell.env = env
        self.env = env

        self.reactants = cell.reactants
        self.products = cell.products

        self.stdE_RTP = self.forward.stdE_RTP - self.reverse.stdE_RTP
        self.n = n # number of moles of electrons transferred per mole of
          # product
        self.H2O = H2O
        #self.sol = reaction.electrolyte(reactants, products)


    def update_quotient(self, getgamma=True, qconc=False, qmolal=False):
        """Calculate the reaction quotient of the redox reaction.

        Note that the forward reaction is the one which appears to be
        dissociating backwards --- its kind of tricky to get your
        head around. We use the method in MT pp 544--546, which I have
        written up in a more mathematically general manner in my notes.

        NOTE: this assumes only the ions taking part have an activity
        other than 1! In most cases this will be valid but not all.
        """

        if self.cell.all_activities == True and not qconc and not qmolal:
            # we have the activities of all reagents, so Q can be calculated
              # using the typical expression
            self.cell.update_quotient(qconc=False, qmolal=False)
              # ^ use the default activity calculator
            self.quotient = self.cell.quotient

        else:
            # we do not have the activities of all our constituents, so
              # calculate the quotient from the mean activity coefficients
              # (find them if needed) and the molalities. This assumes
              # the oxidised reaction becomes its standard state and the
              # reduced reaction leaves its std state.
              # See the theory in the faux-documentation
            # get the salt activities
            fwd_a_salt = self.forward.get_a_salt(getgamma=getgamma)
            rvs_a_salt = self.reverse.get_a_salt(getgamma=getgamma)

            # Find the correct ratio to put into the quotient expression
              # for both the forward:
            fwd_ratio = None
            for r1, mr1 in chain(self.forward.reactants.items(),
              self.forward.products/items()):
                for r2, mr2 in chain(self.cell.reactants.items(),
                  self.cell.products.items()):
                    if r1.name == r2.name: # i.e. how many forward reactions
                      # are needed in the cell reaction!
                        fwd_ratio = mr2
              # and reverse reactions
            rvs_ratio = None
            for r1, mr1 in chain(self.reverse.reactants.items(),
              self.reverse.products.items()):
                for r2, mr2 in chain(self.cell.reactants.items(),
                  self.cell.products.items()):
                    if r1.name == r2.name:
                        rvs_ratio = mr2

            self.quotient = ((rvs_a_salt ** rvs_ratio)
              / (fwd_a_salt ** fwd_ratio))




    def update_E(self, getgamma=True, estimateDifferentials=True):
        """Calculate the electrode potential at the environmental
        temperature and pressure, though the latter is a little
        ambiguous.

        If we have neither, then we can use the reaction quotient,
        which requires our reactants have assigned activities.
        """
        self.forward.update_E(estimate=estimateDifferentials)
        self.reverse.update_E(estimate=estimateDifferentials)


        self.stdE = self.forward.stdE - self.reverse.stdE

        # now correct for nonstandard conditions using activities
        self.update_quotient(getgamma=getgamma)
        self.E = self.stdE - ((self.R*self.env.T/(self.n*self.Fcons))
          * math.log(self.quotient))


    def update_molar_gibbs(self):
        """Update the standard molar gibbs free energy of this reaction.
        """
        if self.E == None:
            raise ValueError("Please first update the electrode potential "
              "with any tabulated data you have before trying to calculate "
              "the free energy!")
        self.molar_gibbs = -self.n * self.Fcons * self.E

    def get_equation(self):
        """Return the overall cell reaction"""
        return self.equation

    def react(self, n):
        """Perform a reaction, consuming unit n moles of reactant.

        Perform the update to the cell reaction. As all of the reactants
        are shared, this should automatically update the fwd and reverse
        reactions too.
        """
        #  this can be looked into in more detial later, if we wish to
          # persue with the redox class. It would be a nice idea to update
          # the gammas after performing the reaction.
        self.cell.react(n)
