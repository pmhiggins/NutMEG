"""
Most recent changes: v2 overhaul 2026

TODO: better handling of phase_ss

@author P M Higgins
@version 0.0.3

"""
import warnings
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
    mol : float
        Total amount of moles present.
    molarity : float
        molarity in mol/L. If gaseous, describes partial pressure in bar.
    molality : float
        molality in mol/kg.
    gamma_molal : float
        Molal activity coefficient
    gamma_molar : float
        Molar activirty coefficient
    activity : float
        Activity
    charge : float
        Charge of reagent. Default 0.
    thermo : bool
        Determines whether or not to compute thermodynamic parameters
        with the host reactor's database.
    phase : str
        identifier for phase of the reagent. Must be one of 'aq','g','s','l'.
    phase_ss : bool
        Identifies whether or not the reagent is in it's standard state.
        Primarily used to know whether to exclude species like H2O from
        quotient calculations.
    """


    def __init__(self, name, locale, amount=None, activity=None,
      gamma_molal=None, gamma_molar = None,
      charge=None, phase='aq', phase_ss=False,
      thermo=True, add_to_locale=True):
        """
        Parameters
        ----------
        name : str
            name of the reagent
        locale : reactor
            Reactor to initialise this reagent in.
        amount : tuple, optional
            Tuple in form (float, str), where the str is an identifier, such as
            'mol', 'molar', 'molal', and the float is the amount value.
        activity : float, optional
            Activity of the species.
        gamma_molal : float, optional
            Molal activity coefficient
        gamma_molar : float, optional
            Molar activirty coefficient
        charge : float, optional
            Charge of reagent. Default 0.
        thermo : bool, optional
            Determines whether or not to compute thermodynamic parameters
            with the host reactor's database.
        phase : str, optional
            identifier for phase of the reagent. Must be one of 'aq','g','s','l'.
        phase_ss : bool, optional
            Identifies whether or not the reagent is in it's standard state.
            Primarily used to know whether to exclude species like H2O from
            quotient calculations.
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


        self.gamma_molal = gamma_molal
        self.gamma_molar = gamma_molar

        self.amount_handler(amount, locale) # sets self.mol, if amount is passed
        self.update_amount(locale, self.mol, activity=activity, zero_warn=False)

        self.charge = charge
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

    def get_amount(self):
        """Return amount in moles"""
        return self.mol

    def get_molality(self):
        """Return molality in mol/kg"""
        return self.molality

    def get_molarity(self):
        """Return molarity in mol/L"""
        return self.molarity

    def get_activity(self):
        """Return activity"""
        return self.activity

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

    def amount_handler(self, amount, locale):
        """
        Parse initialised amount tuple to assign molar amount of species.
        """
        if not amount:
            self.mol = None
            return

        val = amount[0]
        key = amount[1]
        key_lc = amount[1].lower()

        if key_lc =='mole' or key_lc =='mol' or key_lc == 'moles':
            self.mol = val
        elif key == 'm' or key_lc == 'molality' or key_lc == 'molal':
            self.mol = val * locale.kgH2O
        elif key == 'M' or key_lc == 'molarity' or key_lc == 'molar':
            self.mol = val * locale.V_L
        else:
            raise ValueError('Unknown species amount identifier: '+key)


    def update_amount(self, locale,
      mol=None, gamma_molar=None, gamma_molal=None, activity=None,
      zero_warn=True):
        """
        Updates the molar amount of species, including molality and molarity
        if possible. Pass all up-to-date quantities for the best translation.
        """

        if self.name == 'H2O(aq)':
            if activity:
                self.activity = activity
                self.mol = 55.5*locale.kgH2O
            else:
                self.activity = 1.
            return self.mol, self.activity

        # update activity coefficeints if they have been passed.
        if gamma_molar:
            self.gamma_molar = gamma_molar
        if gamma_molal:
            self.gamma_molal = gamma_molal

        if not mol and activity:
            self.activity = activity

            if self.gamma_molal:
                self.molality = activity / self.gamma_molal
                self.mol = self.molality * locale.kgH2O
                if self.gamma_molar:
                    self.molarity = activity / self.gamma_molar
                else:
                    self.molarity = self.mol / locale.V_L

            elif self.gamma_molar:
                self.molarity = activity / self.gamma_molar
                self.mol = self.molarity * locale.V_L
                self.molality = self.mol / locale.kgH2O
            else:
                # neither gammas are known, so functionally assume they are 1.
                # to get an estimated molar quantity.
                self.mol = self.activity * locale.kgH2O

        elif not activity and mol:
            self.mol = mol

            self.molality = self.mol / locale.kgH2O
            self.molarity = self.mol / locale.V_L

            if self.gamma_molal:
                self.activity = self.gamma_molal * self.molality
            elif self. gamma_molar:
                self.activity = self.gamma_molar * self.molality
            else:
                # should we keep this, or refrain from setting activity altogether?
                self.activity = 1. * self.molarity

        elif activity and mol:
            self.mol = mol
            self.activity = activity
            self.molality = self.mol / locale.kgH2O
            self.molarity = self.mol / locale.V_L

        else:
            if zero_warn:
                warnings.warn('no amounts passed to update '+self.name+' with. Setting amounts to zero.')
            self.mol = 1e-16
            self.activity=1e-16

        return self.mol, self.activity





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
