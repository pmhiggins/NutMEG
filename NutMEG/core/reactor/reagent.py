"""
Most recent changes: v2 overhaul 2026

@author P M Higgins
@version 0.0.3

"""

# import NutMEG.util.NutMEGparams as nmp
# from NutMEG.util.loggersetup import loggersetup as logset
# logger = logset.get_logger(__name__, filelevel=nmp.filelevel, printlevel=nmp.printlevel)


class Reagent:
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
    thermo : bool
        Determines whether or not to compute thermodynamic parameters
        with the host reactor's database.
    phase : str
        identifier for phase of the reagent. Must be one of 'aq','g','s','l'.
    phase_ss : bool
        Identifies whether or not the reagent is in it's standard state.
    """


    def __init__(self, name, locale, conc=None, activity=None,
      molal=None, charge=None, gamma=1., phase='aq', phase_ss=False,
      thermo=True, add_to_locale=True):
        """
        Parameters
        ----------
        name : str
            name of the reagent
        conc : float, optional
            molarity in mol/L. If gaseous, conc describes pressure in bar.
            Can be None if molal or activity are known
        gamma : float, optional
            activity coefficient
        activity : float, optional
            Activity. If None, estimate using gamma and conc.
        molal : float, optional
            molality in mol/kg.
            Can be None is conc or activity are known
        charge : float, optional
            Charge of reagent. Default 0.
        thermo : bool, optional
            Determines whether or not to compute thermodynamic parameters
            with the host reactor's database.
        phase : str, optional
            identifier for phase of the reagent. Must be one of 'aq','g','s','l'.
        phase_ss : bool, optional
            Identifies whether or not the reagent is in it's standard state.
        add_to_locale : bool, optional
            Identifies if upon initialisation the reagent should be added to
            the locale's composition. Default True.

        """
        self.name = name

        if phase != 'aq' and phase != 's' and phase != 'g' and phase != 'l':
            raise ValueError("Incorrectly defined phase for reagent "
              + str(name) + ", must be one of 's', 'l', 'g', or 'aq'.")

        # thermodynamic quantities at RTP
        self.std_formation_enthalpy_RTP = None # J/mol
        self.std_formation_entropy_RTP = None # J/mol K
        self.std_formation_gibbs_RTP = None # J/mol
        self.Cp_RTP = None # specific heat capacity

        # thermodynamic quantities at the current evironment
        self.std_formation_gibbs_env = None
        self.std_formation_entropy_env = None
        self.std_formation_enthalpy_env = None
        self.Cp_env = None

        self.thermo = thermo
        self.rkt_twin = None

        # pass Thermo as False to manually update thermochemical parameters.
        # otherwise, set up a reaktoro 'twin'
        if name != 'e-' and self.thermo:
            self.rkt_twin = locale.thermodb.species(self.name)
            self.update_thermo_RTP()
            self.update_thermo(locale)

        elif name == 'e-':
            self.std_formation_enthalpy_RTP = 0.    # J/mol
            self.std_formation_entropy_RTP = 0.    # J/mol K
            self.std_formation_gibbs_RTP = 0.
            self.std_formation_gibbs_env = 0.
            self.std_formation_entropy_env = 0.
            self.std_formation_enthalpy_env = 0.

        self.conc = conc
        self.activity = activity
        self.charge = charge
        self.molal = molal
        self.gamma = gamma
        self.phase = phase
        self.phase_ss = phase_ss

        if add_to_locale:
            locale.add_reagent(self)


    def __str__(self):
        return self.name

    @staticmethod
    def get_phase_str(namestr):
        """Return the phase of a reagent from its name: eg 'aq' from Al(aq)"""
        if namestr[-2]=='q':
            return 'aq'
        elif namestr[-2] == 's' or namestr[-2] == 'l' or namestr[-2] == 'g':
            return namestr[-2]
        else:
            # phase unknown, assume it is aqueous
            return 'aq'


    """

        THERMODYNAMIC CALCULATIONS

    """

    def update_thermo(self, locale):
        """
        Import the reagent's thermal parameters in the current environment.
        """
        if locale.T != 298.15 or locale.P != 101325.0:
            rprops = self.rkt_twin.props(locale.T, 'K', locale.P, 'Pa')
            self.std_formation_gibbs_env = rprops.G0
            self.std_formation_enthalpy_env = rprops.H0
            self.std_formation_entropy_env = rprops.S0
            self.Cp_env = rprops.Cp0
        else:
            # we're in RTP so no need to look up the data again
            self.std_formation_enthalpy_env = self.std_formation_enthalpy_RTP
            self.std_formation_entropy_env = self.std_formation_entropy_RTP
            self.std_formation_gibbs_env = self.std_formation_gibbs_RTP
            self.Cp_env = self.Cp_RTP


    def update_thermo_RTP(self):
        """Import the reagent's thermal parameters at RTP."""
        rprops = self.rkt_twin.props(298.15, 'K', 101325.0, 'Pa')
        self.std_formation_gibbs_RTP = rprops.G0
        self.std_formation_enthalpy_RTP = rprops.H0
        self.std_formation_entropy_RTP = rprops.S0
        self.Cp_RTP = rprops.Cp0


    """

        SETS FOR UPDATING PARAMETERS

    """

    def set_concentration(self, newconc):
        """Update reagent concentration to be newconc in mol/L.

        Parameters
        ----------
        newconc : float
            New molarity to set in mol/L
        """
        self.conc = newconc

    def set_molality(self, newmolal):
        """Update molality in mol/kg solvent.

        Parameters
        ----------
        newmolal : float
            New molality to set.
        """
        self.molal = newmolal

    def set_activitycoefficient(self, newg):
        """Update the activity coefficients.


        Parameters
        ----------
        newg : float
            New activity coefficient to set.
        """
        self.gamma = newg


    def set_phase(self, locale, newphase):
        """
        Change the reagent's phase, and update thermodynamic parameters
        accordingly.


        Parameters
        ---------
        newphase : str
            New phase to set. Must be one of 'aq'. 's', 'g', 'l'
        """
        if phase != 'aq' and phase != 's' and phase != 'g' and phase != 'l':
            raise ValueError("Incorrectly defined phase for reagent "
              + str(name) + ", must be one of 's', 'l', 'g', or 'aq'.")
        self.phase = newphase
        if name != 'e-' and self.thermo:
            self.GetThermoParams(locale)
