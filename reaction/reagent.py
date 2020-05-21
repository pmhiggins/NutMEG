"""

This is the reagent submodule. reagent essentially stores parameters of
individual molecules in a given environment, to be used throughout the
reaction module. Thermodynamic properties are calculated using reaktoro
if available.

Most recent changes: Thermodynamics management December 2019

@author P M Higgins
@version 0.0.3

"""

from NutMEG.environment import environment
from NutMEG.reaction.thermo.reagent_thermo import reagent_thermo

import numpy as np
from os import system
import os.path
import sqlite3

import NutMEG.util.NutMEGparams as nmp
from NutMEG.util.loggersetup import loggersetup as logset
logger = logset.get_logger(__name__, filelevel=nmp.filelevel, printlevel=nmp.printlevel)


class reagent:
    """
    Class for storing and calculating individual reagent properties such
    as concentrations, activities, etc.

    Attributes
    ------------
    name : str
        name of the reagent
    conc : float
        molarity in mol/L. If gaseous, conc describes pressure in bar.
        Can be None if molal or activity are known
    gamma : float
        activity coefficient
    activity : float
        Activity. If None, estimate using gamma and conc.
    molal : float
        molality in mol/kg.
        Can be None is conc or activity are known
    charge : float
        Charge of reagent. Default 0.

    """
    #name = '' #: name of the reagent
    conc = None # mol/l
      # for a gaseous reagent in gaseous reactions, conc describes
      # pressure (in bar).
    gamma = 1. # activity coefficient
    activity = None
    molal = None # molality in mol/kg

    charge = 0

    radius = None #ionic radius used for electrolytes

    phase = None # must be one of 'aq', 'g', 'l', or 's' when initialised
    phase_ss = False # is the reactant in its standard state?

    Cp_RTP = None # specific heat capacity at 298.15 K 100000 Pa
    Cp_env = None
    Cp_T_poly = None # specific heat capacity as a polynomial of temperature
     # type(np.poly1d)

    std_formation_enthalpy_RTP = None # J/mol
    std_formation_entropy_RTP = None # J/mol K
    std_formation_gibbs_RTP = 0. # J/mol

    # thermodynamic quantities at the current evironment
    env = environment(T=298.15, P=101325.0) # use RTP as the default
    std_formation_gibbs_env = None
    std_formation_entropy_env = None
    std_formation_enthalpy_env = None


    # booleans of state
    thermo = True # whether we have the thermodynamic data available


    ####   INITIALISATION METHODS


    def __init__(self, name, env, thermo=True, conc=None, activity=None,
      molal=None, charge=0, gamma=1., radius=None, phase='aq', phase_ss=False,
      Cp_T_poly=None, new=True):
        if new:
            logger.info('Initialising ' + name)
        self.name = name
        self.env = env
        if phase != 'aq' and phase != 's' and phase != 'g' and phase != 'l':
            raise ValueError("Incorrectly defined phase for reagent "
              + str(name) + ", must be one of 's', 'l', 'g', or 'aq'.")
        # pass Thermo as False to update thermochemical parameters yourself
        if name != 'e-' and thermo:
            self.GetThermoParams()
        elif name == 'e-':
            self.std_formation_enthalpy_RTP = 0.    # J/mol
            self.std_formation_entropy_RTP = 0.    # J/mol K
            self.std_formation_gibbs_env = 0.
            self.std_formation_entropy_env = 0.
            self.std_formation_enthalpy_env = 0.
        self.thermo = thermo
        self.conc = conc
        self.activity = activity
        self.charge = charge
        self.molal = molal
        self.gamma = gamma
        self.phase = phase
        self.phase_ss = phase_ss
        self.Cp_T_poly = Cp_T_poly
        self.radius = radius


    def __str__(self):
        return self.name

    @staticmethod
    def get_phase_str(namestr):
        """Get the phase of a reagent from its name: eg 'aq' from Al(aq)"""
        if namestr[-2]=='q':
            return 'aq'
        elif namestr[-2] == 's' or namestr[-2] == 'l' or namestr[-2] == 'g':
            return namestr[-2]
        else:
            # phase unknown, assume it is aqueous
            return 'aq'


    def redefine(self, re):
        """re-initialise as a new or updated reagent"""
        logger.debug('Redefining '+ self.name)
        self.name = re.name
        self.env = re.env
        self.std_formation_gibbs_RTP = re.std_formation_gibbs_RTP
        self.std_formation_enthalpy_RTP = re.std_formation_enthalpy_RTP
        self.std_formation_entropy_RTP = re.std_formation_entropy_RTP
        self.std_formation_gibbs_env = re.std_formation_gibbs_env
        self.std_formation_entropy_env = re.std_formation_entropy_env
        self.std_formation_enthalpy_env = re.std_formation_enthalpy_env
        self.Cp_env = re.Cp_env
        self.Cp_RTP = re.Cp_RTP
        self.thermo = re.thermo
        self.conc = re.conc
        self.activity = re.activity
        self.charge = re.charge
        self.molal = re.molal
        self.gamma = re.gamma
        self.phase = re.phase
        self.phase_ss = re.phase_ss
        self.Cp_T_poly = re.Cp_T_poly
        self.radius = re.radius

        # self = re

    def GetThermoParams(self):
        """Import the reagent's thermal parameters at both RTP and in
        the current environment.
        """
        if self.std_formation_entropy_RTP is None:
            self.import_RTP_params()
        if self.env.T != 298.15 or self.env.P != 101325.0:
            self.import_params_db()
        else:
            # we're in RTP so no need to look up the data again
            self.std_formation_enthalpy_env = self.std_formation_enthalpy_RTP
            self.std_formation_entropy_env = self.std_formation_entropy_RTP
            self.std_formation_gibbs_env = self.std_formation_gibbs_RTP
            self.Cp_env = self.Cp_RTP


    def import_RTP_params(self):
        """Import thermodynamic RTP data from the SQLite database TPdb.
        If the data doesn't exist, calculate it using reaktoro.
        """
        rt = reagent_thermo(self)
        try:
            ThermoData = rt.db_select(T=298.15, P=101325.0)
            self.std_formation_gibbs_RTP = float(ThermoData[0])
            self.std_formation_enthalpy_RTP = float(ThermoData[1])
            self.std_formation_entropy_RTP = float(ThermoData[2])
            self.Cp_RTP = float(ThermoData[3])
        except Exception as e:
            # looks like it isn't in the database, better add it!
            logger.info('Adding ' + self.name
              + ' properties at RTP to the database')

            datasent = rt.thermo_to_db(T=298.15, P=101325.0)

            # now try again if it went in
            if datasent:
                ThermoData = rt.db_select(T=298.15, P=101325.0)
                self.std_formation_gibbs_env = float(ThermoData[0])
                self.std_formation_enthalpy_env = float(ThermoData[1])
                self.std_formation_entropy_env = float(ThermoData[2])
                self.Cp_env = float(ThermoData[3])
            else:
                # this thing shouldn't be using thermo params
                self.thermo = False


    def import_params_db(self):
        """Import thermodynamic data from the SQLite database TPdb.
        If the data doesn't exist, calculate it using reaktoro.
        """
        rt = reagent_thermo(self)
        try:

            ThermoData = rt.db_select()
            self.std_formation_gibbs_env = float(ThermoData[0])
            self.std_formation_enthalpy_env = float(ThermoData[1])
            self.std_formation_entropy_env = float(ThermoData[2])
            self.Cp_env = float(ThermoData[3])
        except Exception as e:
            # looks like it isn't in the database, better add it!
            logger.info('Adding ' + self.name + ' properties at T = '
              + str(round(self.env.T,2)) + ' K and P = ' + str(round(self.env.P,0))
              + ' Pa to the database...')


            datasent = rt.thermo_to_db()

            if datasent:
                # now try again
                ThermoData = rt.db_select()
                self.std_formation_gibbs_env = float(ThermoData[0])
                self.std_formation_enthalpy_env = float(ThermoData[1])
                self.std_formation_entropy_env = float(ThermoData[2])
                self.Cp_env = float(ThermoData[3])
            else:
                # this thing shouldn't be using thermo params
                self.thermo = False




    """

        SETS FOR UPDATING PARAMETERS

    """

    def update_reagent(self): # more to come I'm sure
        """Update the reagent's parameters based on a changing environment.

        For now, it just updates the thermodynamic data.
        """
        if name != 'e-' and thermo:
            self.GetThermoParams()



    def set_concentration(self, newconc):
        """Update reagent concentration in mol/L.
        """
        self.conc = newconc

    def set_molality(self, newmolal):
        """Update molality in mol/kg solvent.

        Parameters
        ---------
        newmolal : float
            New molality to set.
        """
        self.molal = newmolal

    def set_activitycoefficient(self, newg):
        """Update the activity coefficients.
        """
        self.gamma = newg

    def set_phase(self, newphase):
        """Change the reagent's phase, and update thermodynamic parameters
        accordingly.
        """
        if phase != 'aq' and phase != 's' and phase != 'g' and phase != 'l':
            raise ValueError("Incorrectly defined phase for reagent "
              + str(name) + ", must be one of 's', 'l', 'g', or 'aq'.")
        self.phase = newphase
        if name != 'e-' and thermo:
            self.GetThermoParams()
