"""
Most recent changes: v2 overhaul 2026

@author P M Higgins
"""
# import sys
# sys.path.append('../../..')
# from NutMEG.environment import environment
# from NutMEG.reaction.thermo.reaction_thermo import reaction_thermo

import sys
import math
from itertools import chain
import numpy as np
import warnings
from os import system
import os.path
# import sqlite3
from uncertainties import ufloat, umath

gas_const = 8.314472  # J/mol.K

# import NutMEG.util.NutMEGparams as nmp
# from NutMEG.util.loggersetup import loggersetup as logset
# logger = logset.get_logger(__name__, filelevel=nmp.filelevel, printlevel=nmp.printlevel)

class Reaction:

    """Generalised reaction class for unclassified chemical interactions.

    General thermodynamic calculations, and can be used as a basis for child
    modules. Primarily this is for standard conditions only, but for ideal
    fluids discretion for temperature is included, using the reaktoro package.
    If reagents are in the SUPCRT07 database, then free energies can be
    calculated if the reaction proceeds purely thermochemically.

    Attributes
    ----------
    reactants : dict
        Participating reactants in the form {name:molar ratio}
    products : dict
        Participating products in the form {name:molar ratio}
    equation : str
        str of the equation in a readable form.
    frequency_factor : float
        the pre-exponential factor in an arrhenius equation. For rate calculation.
    molar_activation_E : float
        Molar activation energy for an arrhenius equation. For rate calculation.
        Unit J/K mol
    mass_activation_E : float
        Mass activation energy for an arrhenius equation. For rate calculation.
        Unit J/kg
    std_molar_gibbs : float
        Standard Gibbs free energy of reaction in J/mol
    molar_gibbs :
        Molar Gibbs free energy of reaction in J/mol
    mass_gibbs : float
        Gibbs free energy of reaction in J/kg. Sometimes also referred to as
        Gibbs free energy density.
    std_molar_enthalpy :float
        Standard enthalpy of reaction in J/mol
    std_molar_entropy : float
        Standard entropy of reaction in J/mol
    lnK : float
        Natural logarithm of the equilibrium constant. T and P sensitive.
    quotient : float
        Reaction quotient [prod]/[react]
    thermo : bool
        Determines whether or not to compute thermodynamic parameters
        with the host reactor's database.
    """


    def __init__(self, locale, reactants, products,
      frequency_factor=None,
      molar_activation_E=None,
      add_to_locale = True):
        """
        Parameters
        ----------
        locale : reactor
            Host reactor.
        reactants : dict
            Participating reactants in the form {name:molar ratio}
        products : dict
            Participating products in the form {name:molar ratio}
        frequency_factor : float, optional
            the pre-exponential factor in an arrhenius equation. Default None
        molar_activation_E : float
            Molar activation energy for an arrhenius equation. Unit J/K mol.
            Default None.
        add_to_locale : bool, optional
            Identifies if upon initialisation the reagent should be added to
            the locale's composition. Default True.
        """


        self.reactants = reactants
        self.products = products
        self.quotient = None
        self.frequency_factor = frequency_factor
        self.molar_activation_E = molar_activation_E
        self.rate_constant_RTP = None
        self.rate_constant_env = None

        self.all_activities = self.activity_finder(locale)
        self.equation = self.get_equation()
        self.thermo = self.thermo_finder(locale)
        self.rkt_twin = None
        self.lnK = None
        self.stdG = None
        if self.thermo:
            self.rkt_twin = locale.thermodb.reaction(self.equation)
            self.update_thermo(locale)

        if add_to_locale:
            locale.add_reaction(self)


    def __str__(self):
        return self.equation


    def activity_finder(self, locale):
        """Return True if the activities of all of the reagents is known."""
        for r in chain(list(self.products), list(self.reactants)):
            _r = locale.composition[r]
            if _r.activity is None:
                return False # we don't have all the activities
        return True


    def get_equation(self):
        """Return the equation as a string."""
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

        Parameters
        ----------
        r : str
            reactant identifier (e.g. 'H2O(l)')
        mr : float
            molar ratio in the reaction
        """
        result = ""
        if mr != 1:
            result += str(mr) + "*" + str(r) + " "
        else:
            result += str(r) + " "
        return result



    def thermo_finder(self, locale):
        """
        Return True if thermodynamic data is available and calculable
        for all reagents.
        """
        for re in chain(list(self.reactants), list(self.products)):
            _re = locale.composition[re]
            if _re.name !='e-' and not _re.thermo:
                return False # At least one does not have the data
        # we must have everything
        return True



    ####### SETS FOR UPDATING PARAMETERS  ######


    def set_reactants(self, re):
        """Redefine dictionary of reactants as re."""
        self.reactants = re


    def set_products(self, pr):
        """Redefine dictionary of products as pr."""
        self.products = pr



    ########  GENERIC CALCULATIONS: RATES


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
              math.exp(-self.molar_activation_E/(gas_const*self.env.T)))




    ####### GENERIC CALCULATIONS: THERMODYNAMICS


    def quotient_calculator(self, locale, attr):
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
            _p = locale.composition[p]
            if _p.phase_ss == False:
                try:
                    A = float(getattr(_p, attr))
                    if A != 0.:
                        a = float(mr)
                        multiplier = multiplier * math.pow(A, a)
                except:
                    if attr=='activity':
                        A=_p.activity.n
                        # print(p, A)
                    if A != 0.:
                        a = float(mr)
                        multiplier = multiplier * umath.pow(A, a)

        for r, mr in self.reactants.items():
            _r = locale.composition[r]
            if _r.phase_ss == False:
                try:
                    A = float(getattr(_r, attr))
                    if A != 0.:
                        a = float(mr)
                        multiplier = multiplier / math.pow(A, a)
                except:
                    if attr=='activity':
                        A=_r.activity.n
                        # print(r, A)
                    if A != 0.:
                        a = float(mr)
                        multiplier = multiplier * umath.pow(A, a)
        return multiplier




    def update_quotient(self, locale, qconc=False, qmolal=False):
        """Update the reaction quotient for this reaction.

        Parameters
        ----------
        qconc : bool, optional
            If True, calculate using molarity (default is False).
        qmolal : bool, optional
            If True, calculate using molality (default is False).

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
            multiplier = (self.quotient_calculator(locale, "conc")
              * self.quotient_calculator(locale, "gamma"))
        elif qmolal==True:
            # Calculate using molality
            multiplier = (self.quotient_calculator(locale, "molal")
              * self.quotient_calculator(locale, "gamma"))
        else:
            # The default is to use ativities
            multiplier = self.quotient_calculator(locale, "activity")
        self.quotient = multiplier




    def update_std_molar_enthalpy_of_reaction(self, locale):
        """Update the standard molar enthalpy of reaction from the
        enthalpies of formation of the reagents in the current
        environment.
        """

        HoR = 0.
        for p, mr in self.products.items():
            _p = locale.composition[p]
            HoR += (_p.std_formation_enthalpy_env * mr)
        for r, mr in self.reactants.items():
            _r = locale.composition[r]
            HoR -= (_r.std_formation_enthalpy_env * mr)
        self.std_molar_enthalpy = HoR



    def update_std_molar_entropy_of_reaction(self, locale):
        """Update the standard molar entropy of reaction from the
        entropies of formation of the reagents in the current
        environment.
        """
        SoR = 0.
        for p, mr in self.products.items():
            _p = locale.composition[p]
            SoR += (_p.std_formation_entropy_env * mr)
        for r, mr in self.reactants.items():
            _r = locale.composition[r]
            SoR -= (_r.std_formation_entropy_env * mr)
        self.std_molar_entropy = SoR




    def update_std_molar_gibbs_from_reagents(self, locale):
        """
        Update the standard molar gibbs free energy of reaction
        using the Gibbs free energy of formation of the reagents in
        the current environment, if available.

        It is preferable to do this using reaktoro (update_thermo())
        if the thermodynamic data is available in the standard databases.
        """
        GoR = 0.
        for p, mr in self.products.items():
            _p = locale.composition[p]
            GoR += (_p.std_formation_gibbs_env * mr)
        for r, mr in self.reactants.items():
            _r = locale.composition[r]
            GoR -= (_r.std_formation_gibbs_env * mr)

        self.std_molar_gibbs = GoR




    def update_molar_gibbs_from_quotient(self, locale,
      Q_qconc=False,
      Q_qmolal=False):
        """
        Update Gibbs free energy of reaction at temperature T,
        using the expression:
        :math:`\Delta G_{T} = \Delta G_{T}^{0} + RT\ln{Q}`

        #TODO: improve error handling.
        """

        # update Q and proceed
        self.update_quotient(locale, qconc=Q_qconc, qmolal=Q_qmolal)

        try:
            self.molar_gibbs = (self.std_molar_gibbs
              + (gas_const * locale.T * math.log(self.quotient)))
        except:
            self.molar_gibbs = None


    def update_mass_gibbs(self, locale):
        """
        Return an approximation of the energy density [J/kg H2O] for the
        reaction

        Calculates the smallest energy yield from 'using up' the
        reagents. In reality, the free energy would change as the concentration
        decreases, so this is only a measure of the energy density available
        for this reaction at this moment in time.
        """
        ED = []
        for r, mr in self.reactants.items():
            _r = locale.composition[r]
            if _r.name != 'H2O(aq)' and _r.name != 'H+' and _r.name != 'OH-':
                ED.append(_r.molal*-self.molar_gibbs/mr)
        self.mass_gibbs = min(ED)


    def react(self, n, locale):
        """Perform a reaction, consuming unit n moles of reactants.

        Parameters
        ----------
        n : float
            Total number of moles to react throughout the reactor
        locale : reactor
            The reactor containing this reaction's reagents.

        Notes
        -----
        n is the number of moles of the reaction occuring, so if both
        reactants had a molar ratio of 4, and n=1 was passed, 4 moles
        of each reactant would be consumed.

        # TODO: would it be better to focus on total moles, so V and kgH2O
        are less prominent?
        """
        for r, mr in self.reactants.items():
            _r = locale.composition[r]
            # Find total number of moles in system, then remove the amount
            # that has been reacted away or formed.
            _r.conc = (((_r.conc*1000.0*locale.V)
              - (mr*n))/(1000.0*locale.V))
            _r.molal = (((_r.molal*locale.kgH2O)
              - (mr*n))/(locale.kgH2O))
            if _r.conc<0:
                _r.conc=0
            if _r.molal<0:
                _r.molal=0
            if _r.name != 'H2O(l)':
                _r.activity = _r.conc * _r.gamma
        for p, mr in self.products.items():
            _p = locale.composition[p]
            # Find total number of moles in system, then remove the amount
            # that has been reacted away or formed.
            _p.conc = (((_p.conc*1000.0*locale.V)
              + (mr*n))/(1000.0*locale.V))
            _p.molal = (((_p.molal*locale.kgH2O)
              + (mr*n))/(locale.kgH2O))
            if _p.name != 'H2O(l)':
                _p.activity = _p.conc * _p.gamma



    ####### UPDATING USING REAKTORO

    ###### Use the reaktoro package to perform thermodynamic calculations


    def update_thermo_reagents(self, locale):
        """Update the energetic parameters of the reagents using reaktoro.
        """
        for r in chain(list(self.reactants), list(self.products)):
            _r = locale.composition[r]
            if r.name != 'e-' and r.name != 'H+':
                r.update_thermo(locale)

    def update_thermo(self, locale):
        """Get useful energetic parameters (standard molar Gibbs and lnK)
        for the current state of this reaction.
        """

        self.rkt_twin = locale.thermodb.reaction(self.equation)
        rprops = self.rkt_twin.props(locale.T, 'K', locale.P, 'Pa')

        # update reaction parameters
        self.std_molar_gibbs = rprops.dG0
        self.lnK = rprops.lgK * math.log(10)
