"""

This is the reaction module, which contains the reaction class.
First initialise reagents and then initialise a basic reaction ---
if further specialties are required such as redox,
import one of those submodules (if included).

Most recent changes: Thermodynamics management December 2019

@author P M Higgins
@version 0.1.0

"""

from NutMEG.environment import environment
from NutMEG.reaction.thermo.reaction_thermo import reaction_thermo

import sys
import math
from itertools import chain
import numpy as np
import warnings
from os import system
import os.path
import sqlite3

R = 8.314472  # J/mol.K

import NutMEG.util.NutMEGparams as nmp
from NutMEG.util.loggersetup import loggersetup as logset
logger = logset.get_logger(__name__, filelevel=nmp.filelevel, printlevel=nmp.printlevel)

class reaction:

    """Generalised reaction class for unclassified chemical interactions.

    General thermodynamic calculations, and can be used as a basis for child
    modules. Primarily this is for standard conditions only, but for ideal
    fluids discretion for temperature is included, using the reaktoro package.
    If our reagents are in the SUPCRT07 database, then free energies can be
    calculated if the reaction proceeds purely thermochemically.

    """

    reactants = None  # dict of reagents in the form {reactant:molar ratio}
    products = None  # as above, for products
    equation = ""  #  str of the equation in a readable form.

    # rate parameters
    rate_constant_RTP = None  # units will vary, at RTP
    rate_constant_env = None # value in the local environment
    frequency_factor = None  # the A term in an arrhenius equation
    molar_activation_E = None  # activation energy in J/K mol
    mass_activation_E = None  # activation energy in J/kg

    # energetic parameters
    std_molar_gibbs = None	 # Standard Gibbs free energy of reaction in J/mol
    molar_gibbs = None  # molar Gibbs free energy of reaction in J/mol
    mass_gibbs = None  # Gibbs free energy of reaction in J/kg
    std_molar_enthalpy = None  # Standard enthalpy of reaction in J/mol
    std_molar_entropy = None  # Standard entropy of reaction in J/mol
    lnK = None  # Natural logarithm of the equilibrium constant
    quotient = None  # Reaction quotient [prod]/[react]

    # booleans for the state of the reaction
    equilibrium = False	 # Whether the reaction is at equilibrium or not.
    all_activities = False  # Whether we have activities for all the reagents.

    env = None # The environment.

    """

    INITIALISATION

    """

    def __init__(self, reactants, products, env, frequency_factor=0.,
    molar_activation_E=0., equilibrium=False):
        self.reactants = reactants
        self.products = products
        self.env = env
        self.frequency_factor = frequency_factor
        self.molar_activation_E = molar_activation_E
        self.equilibrium=equilibrium
        self.all_activities = self.activity_finder()
        self.equation = self.get_equation()
        self.thermo = self.thermo_finder()

    # def __str__(self):
    #     return self.equation

    def activity_finder(self):
        """Return True if we have the activities for all of the reagents."""
        for r in chain(list(self.products), list(self.reactants)):
            if r.phase == 's' or r.phase == 'l' or r.name == 'e-':
                # set that activity to 1 if we haven't done so already
                if r.activity == None:
                    r.activity = 1.
                continue
            elif r.activity is None:
                # we clearly don't have all the activities
                return False

        return True


    def get_equation(self):
        """Get the equation as a string."""
        reactionstring = ""

        for r, mr in self.reactants.items():
            reactionstring += self.equationbuilder(r, mr) + "+ "
        else:
            # delete the last plus and turn it into an equals sign.
            # not very elegant but its a simple way to do this.
            reactionstring = reactionstring[:-2] + '= '
        # Now do the same for the products
        for p, mr in self.products.items():
            reactionstring += self.equationbuilder(p, mr) + "+ "
        else:
            reactionstring = reactionstring[:-3]
        return reactionstring


    def equationbuilder(self, r, mr):
        """Return reagent as it would appear in an equation.

        @param r reagent to include
        @param mr its molar ratio"""
        result = ""
        if mr != 1:
            result += str(mr) + "*" + str(r.name) + " "
        else:
            result += str(r.name) + " "
        return result


    def reagents_name(self):
        """Get the reactants and products dictionaries with the names
        as keys rather than the reagent object. """
        reac, prod = {}, {}
        for k, v in self.reactants.items():
            reac[k.name] = v
        for k, v in self.products.items():
            prod[k.name] = v
        return reac, prod


    def thermo_finder(self):
        """Return True if data is available for all reagents in the
        SUPCRT07 Database.
        """
        for re in chain(list(self.reactants), list(self.products)):
            if re.name !='e-' and not re.thermo:
                return False # At least one does not have the data
        return True # We have everything!


    """

    SETS FOR UPDATING PARAMETERS

    """


    def update_reagents(self):
        """Update the parameters of the reagents e.g. when the environ-
        ment has changed, new conc. etc.
        """
        for r in chain(list(self.reactants), list(self.products)):
            r.update_reagent()


    def set_reactants(self, re):
        """Redefine dictionary of reactants as re."""
        self.reactants = re


    def set_products(self, pr):
        """Redefine dictionary of products as pr."""
        self.products = pr



    """

    GENERIC CALCULATIONS: RATES

    """

    def bool_rate_constants(self):
        """Returns bool showing if we have any rate constants."""
        if self.rate_constant_env == None and self.rate_constant_RTP== None:
            return False
        else:
            # We have at least one
            return True

    def calculate_rate(self):
        """Update the rate constant using the arrhenius equation."""

        if self.frequency_factor == None or self.molar_activation_E == None:
            raise ValueError("You have not initialised the frequency_factor "+\
              "and/or the molar_activation_E for your reaction. An "+\
              "Arrhenius calculation cannot be performed.")
        else:
            self.rate_constant_env = (self.frequency_factor *
              math.exp(-self.molar_activation_E/(R*self.env.T)))


    """

    GENERIC CALCULATIONS: THERMODYNAMICS

    """


    def quotient_calculator(self, attr):
        """Return the reaction quotient based on the reactant attribute
        passed.

        Parameters
        ----------
        attr : str
            activity-like attribute of reagent

        Notes
        ----------
        Recognised attrs include: "conc", "molal", "activity".
        """
        multiplier =1.
        A = 1.
        a = 1.
        for p, mr in self.products.items():
            if p.phase_ss == False:
                A = float(getattr(p, attr))
                if A != 0.:
                    a = float(mr)
                    multiplier = multiplier * math.pow(A, a)
        for r, mr in self.reactants.items():
            if r.phase_ss == False:
                A = float(getattr(r, attr))
                if A != 0.:
                    a = float(mr)
                    multiplier = multiplier / math.pow(A, a)
        return multiplier




    def update_quotient(self, qconc=False, qmolal=False):
        """Update the reaction quotient for this reaction.

        Parameters
        ----------
        qconc : bool, optional
            If True, calculate using molarity (default is False).
        qmolal : bool, optional
            If True, calculated using molality (default is False).

        Notes
        ----------
        The default is to use activities, but if needed the optional
        arguments may be switched for using concentrations or molalities
        with activity coefficients.  Valid for gaseous and aqueous
        reagents as conc is defined as equivalent to gas pressure
        in the reagent object.
        Molarity takes precedence over molality.
        """

        multiplier = 1.

        if qconc==True:
            # Calculate using concentrations
            multiplier = (self.quotient_calculator("conc")
              * self.quotient_calculator("gamma"))
        elif qmolal==True:
            # Calculate using molality
            multiplier = (self.quotient_calculator("molal")
              * self.quotient_calculator("gamma"))
        else:
            # The default is to use ativities
            multiplier = self.quotient_calculator("activity")
        self.quotient = multiplier


    def update_std_molar_gibbs_from_quotient(self,
      Q_qconc=False, Q_qmolal=False):
        """Calculate the standard molar gibbs free energy of reaction.

        Notes
        ------
        Uses the expression:
        :math:`\Delta G_{T}^{0} = -RT\ln{K}`

        This can only be done at equilibrium.
        """
        if self.equilibrium == False:
            raise ValueError("This reaction is not at equilibrium, "
            "we cannot use the expression \Delta G0 = -RTln K")
        # Update the equilibrium quotient
        self.update_quotient(qconc=Q_qconc, qmolal=Q_qmolal)
        self.lnK = math.log(self.quotient)
        self.std_molar_gibbs = -(R)*self.env.T*self.lnK


    def update_std_molar_enthalpy_of_reaction(self):
        """Update the standard molar enthalpy of reaction from the
        enthalpies of formation of the reagents in the current
        environment.
        """

        HoR = 0.
        for p, mr in self.products.items():
            HoR += (p.std_formation_enthalpy_env * mr)
        for r, mr in self.reactants.items():
            HoR -= (r.std_formation_enthalpy_env * mr)
        self.std_molar_enthalpy = HoR



    def update_std_molar_entropy_of_reaction(self):
        """Update the standard molar entropy of reaction from the
        entropies of formation of the reagents in the current
        environment.
        """
        SoR = 0.
        for p, mr in self.products.items():
            SoR += (p.std_formation_entropy_env * mr)
        for r, mr in self.reactants.items():
            SoR -= (r.std_formation_entropy_env * mr)
        self.std_molar_entropy = SoR




    def update_std_molar_gibbs_G(self):
        """Update the standard molar gibbs free energy of reaction
        using the Gibbs free energy of formation of the reagents in
        the current environment, if available.

        It is preferable to do this using reaktoro, (rto_current_env)
        if the thermodynamic data is available in the standard databases.
        """
        GoR = 0.
        for p, mr in self.products.items():
            GoR += (p.std_formation_gibbs_env * mr)
            print(p.name, mr, p.std_formation_gibbs_env)
        for r, mr in self.reactants.items():
            GoR -= (r.std_formation_gibbs_env * mr)
            print(r.name, mr, p.std_formation_gibbs_env)
        self.std_molar_gibbs = GoR



    def update_std_molar_gibbs_HS(self):
        """Update the standard molar gibbs free energy of reaction using
        the standard enthalpies and entropies of formation of the reagents,
        if available.

        H has a weak dependence on T, so large extremes will be inresaingly
        poorly represented. S also has a dependence on T, changing with
        heat capacity, latent heat of fusion etc.

        If you do not have the non-RTP values of H or S --- or they are not,
        avaliable in the SUPCRT07 database --- consider using the
        Tdep reaction type (special/Tdep)
        """
        # failsafe in case we have forgotten to calculate the enthalpies
          # and entropies
        if self.std_molar_enthalpy == None:
            self.update_std_molar_enthalpy_of_reaction()
        if self.std_molar_entropy == None:
            self.update_std_molar_entropy_of_reaction()
        GoR = self.std_molar_enthalpy - (self.env.T*self.std_molar_entropy)
        self.std_molar_gibbs = GoR


    def update_molar_gibbs_from_quotient(self, Q_qconc=False,
      Q_qmolal=False, updatestdGibbs=True):
        """Update Gibbs free energy of reaction at temperature T,
        from \DeltaG_T = \DeltaG_T^0 + RTlnQ

        Should demonstrate DeltaG=0 at equilibrium. Again, be cautious
        of \DeltaG0's dependence on T.
        """

        # update Q and proceed
        self.update_quotient(qconc=Q_qconc, qmolal=Q_qmolal)
        if updatestdGibbs: # update the standard Gibbs free energy
            if not self.thermo:
                self.update_std_molar_gibbs_HS()
            else:
                try:
                    self.rto_current_env()
                except:
                    # insufficient info on reactants for reaktoro
                    self.update_std_molar_gibbs_G()
        try:
            self.molar_gibbs = (self.std_molar_gibbs
              + (R * self.env.T * math.log(self.quotient)))
        except:
            self.molar_gibbs = 0

    def react(self, n):
        """Perform a reaction, consuming unit n moles of reactant.

        n is the number of moles of the reaction occuring, so if both
        reactants had a molar ratio of 4, and n=1 was passed, 4 moles
        of each reactant would be consumed.
        """
        for r, mr in self.reactants.items():
            # Find total number of moles in system, then remove the amount
            # that has been reacted away or formed.
            r.conc = (((r.conc*1000.0*self.env.V)
              - (mr*n))/(1000.0*self.env.V))
            if r.conc<0:
                r.conc=0
            if r.name != 'H2O(l)':
                r.activity = r.conc * r.gamma
        for p, mr in self.products.items():
            # Find total number of moles in system, then remove the amount
            # that has been reacted away or formed.
            p.conc = (((p.conc*1000.0*self.env.V)
              + (mr*n))/(1000.0*self.env.V))
            if p.name != 'H2O(l)':
                p.activity = p.conc * p.gamma



    """

    UPDATING USING REAKTORO

    Use the reaktoro package to perform thermodynamic calculations

    """

    def rto_reagents(self):
        """Update the energetic parameters of the reagents using reaktoro.
        """
        for r in chain(list(self.reactants), list(self.products)):
            if r.name != 'e-':
                r.import_params_db()

    def rto_current_env(self):
        """Get useful energetic parameters (standard molar Gibbs and lnK)
        for the current state of this reaction.
        """

        rt = reaction_thermo(self)
        stdG, lnK = rt.get_stdG_lnK()

        # update reaction parameters
        self.std_molar_gibbs = float(stdG)
        self.lnK = float(lnK)
