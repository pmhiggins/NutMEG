"""
Most recent changes: v2 overhaul 2026.

@author P M Higgins
"""
# import sys
# sys.path.append('../../..')

from .reaction import reaction as rxn
from .reagent import reagent as rgt

import warnings
from itertools import chain
from copy import copy
# import sqlite3
import sys, os, ast
# from datetime import date
# from .reactor_dbhelper import rdb_helper

# import NutMEG.util.NutMEGparams as nmp
# from NutMEG.util.loggersetup import loggersetup as logset
# logger = logset.get_logger(__name__, filelevel=nmp.filelevel, printlevel=nmp.printlevel)

class reactor:
    """
    Class for storing physicochemical backdrop of a chemcial environment.
    Contains reagents and reactions, and is able to perform them.

    Attributes
    ----------
    T : float
        Reactor temperature in K
    P : float
        Reactor pressure in Pa
    V : float
        Reactor volume in m^3
    kgH2O : float
        kg of H2O in the reactor
    reactionlist : dict
        reactions which may be performed in the format {'reaction equation' :
        {reaction}}.
    composition : dict
        reagents which can be found in the reactor in the format
        {'reagent name' : reagent}
    pH : float
        pH of the reactor, important for some interactions. Default is None
    composition_inputs : dict, kwarg
        If there is a net flow of reagents in/out of the reactor, add them
        to this dictionary with their name as the key, and rate in M/s as the
        value.
    thermodb : reaktoro database
        Thermodynamic database to be used for all thermodynamic calculations.
    thermodb_name : str
        Name identifier for thermodb.
    """



    def __init__(self, T=298.15, P=101325.0, thermodb='supcrt07-organics', dbtype=None,
      species = [],
      pH = None,
      V = 0.001,
      kgH2O = 1.,
      **kwargs):
        """
        Parameters
        ----------
        T : float, optional
            Reactor temperature in K, Default 298.15
        P : float, optional
            Reactor pressure in Pa, default 101325.0
        V : float, optional
            Reactor volume in m^3, default 0.001
        pH : float, optional
            pH of the reactor. Default is None
        thermodb : 'str', optional
            Name identifier for thermodynamic database to be used, or filename.
        dbtype : str or reaktoro database
            if a string: 'supcrt' or 'pitzer', to denote the type of database
            NutMEG supports. You may also pass a reaktoro database you
            initialise yourself (e.g. NASA) but it may not be supported everywhere.
        composition_inputs : dict, kwarg
            If there is a net flow of reagents in/out of the reactor, add them
            to this dictionary with their name as the key, and rate in M/s as the
            value.
        """

        self.T = T
        self.P = P

        self.thermodb = None
        self.thermodb_name = None
        if isinstance(thermodb, str):
            self.thermodb_name = thermodb
            self.thermodb = self.get_rkt_database(dbtype)


        self.composition = {}
        for s in species:
            if isinstance(s, str):
                self.composition[s] = rgt.from_rkt(self, s)
            else:
                raise ValueError('unknown species object initialised with reactor')

        self.V = V
        self.pH = pH
        self.reactionlist = {}

        self.composition_inputs = kwargs.pop('composition_inputs', {})



    @classmethod
    def from_rto_state(cls, state, props, kgH2O=None):
        """
        Build a reactor object from a reaktoro ChemicalState and its
        ChemicalProps attribute. This method populates a reactor with
        the aqueous and gaseous species present in the reaktoro ChemicalState.

        Parameters
        ----------
        state : reaktoro.ChemicalState
            The source reaktoro chemical state
        props : reaktoro.ChemicalProps
            The source reaktoro chemical props
        kgH2O : float, optional
            To compute molality, pass the kg of liquid H2O present in the
            chemical state. Default None = molalite will not be computed.
        """

        R = cls(name)

        # set reactor physcial env to match the rto state.
        R.T, R.P = float(state.temperature()), float(state.pressure())

        # The reaktoro ChemicalSystem.
        # This manages the phases, which in turn manage the species.
        system = state.system()
        R.thermodb = system.database()

        rto_nm_phases = {'AqueousPhase': 'aq',
          'GaseousPhase': 'g'}
        # mineralphases go by the name of the mineral in rkt,
        # so we need to rethink how to include 's'

        # populate the reactor with reagents read from the reaktoro system
        for phase in system.phases():
            for species in phase.species():
                _n = species.name()
                molal = None
                if kgH2O:
                    molal = float(props.speciesAmount(_n)/kgH2O)

                _s = rgt(_n,
                  R, phase=rto_nm_phases[phase.name()],
                  activity=float(props.speciesActivity(_n)),
                  gamma=float(props.speciesActivityCoefficient(_n)),
                  conc=float(props.speciesConcentration(_n)),
                  molal=molal,
                  charge=float(species.charge())
                )

        return R


    def update_from_rto_state(self, state, props, kgH2O=None):
        """
        Update species in the reactor based on a reaktoro ChemicalState and its
        ChemicalProps attribute. This method only updates existing species
        and won't add new ones. If the reaktoro temperature or pressure has
        changed, so too will the thermodynamic properties of reagents and
        reactions in the reactor.

        Parameters
        ----------
        state : reaktoro.ChemicalState
            The source reaktoro chemical state
        props : reaktoro.ChemicalProps
            The source reaktoro chemical props
        kgH2O : float, optional
            To compute molality, pass the kg of liquid H2O present in the
            chemical state. Default None = molalite will not be computed.
        """

        if self.T != float(state.temperature()) or self.P != float(state.pressure()):
            self.change_TP(T=float(state.temperature()), P=float(state.pressure()))

        # The reaktoro ChemicalSystem.
        # This manages the phases, which in turn manage the species.
        system = state.system()

        rto_nm_phases = {'AqueousPhase': 'aq',
          'GaseousPhase': 'g'}
        # mineralphases go by the name of the mineral in rkt,
        # so we need to rethink how to include 's'

        # populate the reactor with reagents read from the reaktoro system
        for phase in system.phases():
            for species in phase.species():
                _n = species.name()
                if _n in self.composition:
                    self.composition[_n].activity = float(props.speciesActivity(_n))
                    self.composition[_n].gamma = float(props.speciesActivityCoefficient(_n))
                    self.composition[_n].conc = float(props.speciesConcentration(_n))
                    if kgH2O:
                        self.composition[_n].molal = float(props.speciesAmount(_n)/kgH2O)

                    if TPchange:
                        self.change_TP(T=T, P=P)




    def get_rkt_database(self, dbtype=None):
        """
        Return a reaktoro database for extracting thermodynamic parameters.

        Parameters
        ----------
        dbtype : str or reaktoro database
            if a string: 'supcrt' or 'pitzer', to denote the type of database
            NutMEG supports. You may also pass a reaktoro database you
            initialise yourself (e.g. NASA) but it may not be supported everywhere.

        Notes
        -----
        This function also houses lazy imports of the reaktoro database classes
        so only the required one is imported.
        """

        if isinstance(dbtype, str) or not dbtype:
            if self.thermodb_name.startswith('supcrt') or dbtype == 'supcrt':
                from reaktoro import SupcrtDatabase
                return SupcrtDatabase(self.thermodb_name)
            elif dbtype.lower() == 'phreeqc':
                from reaktoro import PhreeqcDatabase
                return PhreeqcDatabase(self.thermodb_name)
            else:
                warnings.warn('Unable to set reaktoro database')
                self.thermodb_name = None
                return None
        else:
            try:
                return dbtype(self.thermodb_name)
            except:
                warnings.warn('Unable to set reaktoro database')
                self.thermodb_name = None
                return None




    def print_reactions(self):
        """ print the list of reaction equations.
        """
        print(self.reactionlist.keys())

    def print_composition(self):
        """print a list of the composition"""
        print(self.composition.keys())



    def add_reaction(self, rxxn):
        """
        Add a new reaction to reactionlist

        Parameters
        ----------
        rxxn : reaction
            Reaction object to be added
        """
        if rxxn.equation in list(self.reactionlist.keys()):
            raise ValueError('Reaction with this eq already present in this reactor!')
        # add the reaction into reationlist
        self.reactionlist[rxxn.equation] = rxxn


    def add_reagent(self, rct):
        """
        Add a new reagent to the reactor composition.

        Parameters
        ----------
        rct : reagent
            Reagent object to be added
        """
        if rct.name in list(self.composition.keys()):
            raise ValueError('New reagent tried to overwrite one in this reactor!')
        # otherwise it's brand new, add away.
        self.composition[rct.name] = rct


    def update_composition(self, t):
        """ If there are inflows into the composition, make the changes there
        would be in time t [s]"""
        for c, r in self.composition_inputs.items():
            self.composition[c].activity += r*t
            self.composition[c].conc += r*t


    def take_step(self, t):
        """
        perform time-sensitive updates to the reactor's attributes with time
        step t. The basic case
        is to update the composition by any fixed rate changes in
        reactor.composition_inputs.

        TODO: consider an option to progress reactions by their rate (if present)
        """
        self.update_composition(t)


    # def unify_reaction(self, rxxn, overwrite=False):
    #     """Add the reaction and its reagents to the reactor, ensuring there is
    #     only one of each reagent type in the reactor.
    #
    #     Returns the reaction after unification
    #
    #     Parameters
    #     ----------
    #     rxxn : NutMEG.reaction.reaction like
    #         reaction to unify
    #     overwrite : bool
    #         if True, overwrite the current composition with the activities of
    #         the reagents in rxxn. Default is False.
    #
    #     Notes
    #     -----
    #     The reaction passed will end up pointing to the composition of the
    #     reactor. Pass overwrite as True to update the values in the composition,
    #     False to use the values we already have. If overwrite is False, and the
    #     reaction has a reagent which is not in this reactor's composition, it
    #     will be added.
    #     """
    #
    #     if not overwrite:
    #         # we want the reaction passed to be reset to be the reaction
    #         # in reactionlist, if it exists in there.
    #         # if not, we need to add it, and unify the reagents with this
    #         # reactor's composition.
    #         try:
    #             # see if the reaction is already saved in the reactor
    #             rxxn = self.reactionlist[rxxn.name][type(rxxn)]
    #             # if it is, return that reaction
    #             return rxxn
    #         except:
    #             logger.info('Reaction to be unified not found in '+self.name)
    #             # return rxxn
    #
    #     if overwrite:
    #         # overwrite the entry in reactionlist with the passed reacrion,
    #         # compositions and all.
    #         self.reactionlist[rxxn.equation][type(rxxn)] = rxxn
    #
    #     # redefine the reagents
    #     # according to the composition.
    #     for rrxn in list(rxxn.reactants.keys()):
    #         inlist = False
    #         for c_name, c_rxt in self.composition.copy().items():
    #
    #             if rrxn.name == c_name:
    #                 inlist=True
    #                 # this reagent is in both the passed reaction and
    #                 # the composition.
    #                 if overwrite:
    #                     # update the composition entry with ragent data from
    #                     # the reaction
    #                     self.composition[c_name].redefine(rrxn)
    #                 # reset the reagent to be the same object as the one
    #                 # in the composition
    #                 rxxn.reactants[c_rxt] = rxxn.reactants.pop(rrxn)
    #         if not inlist:
    #             # it wasn't found in the composition, add it.
    #             # logger.info('Adding '+rrxn.name+' to '+self.name+\
    #             #   "'s composition.'")
    #             self.composition[rrxn.name] = rrxn
    #
    #     for rrxn in list(rxxn.products.keys()):
    #         inlist = False
    #         for c_name, c_rxt in self.composition.copy().items():
    #
    #             if rrxn.name == c_name:
    #                 inlist=True
    #                 # this reagent is in both the passed reaction and
    #                 # the composition.
    #                 if overwrite:
    #                     # update the composition entry with ragent data from
    #                     # the reaction
    #                     self.composition[c_name].redefine(rrxn)
    #                 # reset the reagent to be the same object as the one
    #                 # in the composition
    #                 rxxn.products[c_rxt] = rxxn.products.pop(rrxn)
    #
    #         if not inlist:
    #             # it wasn't found in the composition, add it.
    #             # logger.info('Adding '+rrxn.name+' to '+self.name+\
    #             #   "'s composition.'")
    #             self.composition[rrxn.name] = rrxn
    #     # return rxxn



    def perform_reaction(self, re_eq, n):
        """Perform the reaction n molar times.

        Parameters
        ----------
        re_eq : str
            reaction equation as it appears in reactionlist. If this equation is
            not in reactionlist an error is raised.
        n : float
            number of moles of reaction to perform.
        """
        self.reactionlist[re_eq].react(n)


    def change_TP(self, T=None, P=None, update_thermo=True):
        """
        Update reactor temperature and pressure. If a known thermodynamic
        database is available, standard thermodynamic properties of all hosted
        reagents and reactions will be updated.

        Parameters
        ----------
        T : float
            New temperature in Kelvin
        P : float
            New pressure in Pa
        update_thermo : bool
            Whether to update standard thermodynamic properties of all
            hosted reagents and reactions.
        """
        if T == self.T and P == self.P:
            return
        if T and T != self.T:
            self.T = float(T)
        if P and P != self.P:
            self.P = float(P)
        if update_thermo:
            for k,v in self.composition.items():
                v.update_thermo(self)
            for k,v in self.reactionlist.items():
                v.update_thermo(self)


    def contains_reagent(self, rname):
        """ Return true if a named reagent is present in composition. """
        if rname in self.composition:
            return True
        else:
            return False


    # def Comp_to_db(self, pH, CompID=None, dbpath=nmp.std_dbpath):
    #     """ Send the composition data to the database"""
    #     db = sqlite3.connect(dbpath)
    #     cursor = db.cursor()
    #
    #     compdictstr = self.getconcs()
    #     print(compdictstr)
    #
    #
    #     # fill in the methanogen database
    #     try:
    #         if CompID is None:
    #             #use standard generated compID
    #             cursor.execute(' INSERT INTO composition(CompID, Dict_Pop, pH) VALUES(?,?,?)', ('Tester', compdictstr, pH))
    #             cursor.execute('SELECT rowid FROM Composition WHERE CompID = ?', ('Tester',))
    #             entryno = cursor.fetchone()[0]
    #             CompID = str(entryno)+date.today().strftime("_%d%m%y")
    #             cursor.execute('UPDATE Composition SET CompID = ? WHERE CompID = ?', (CompID, 'Tester'))
    #             self.CompID = CompID
    #         else:
    #             cursor.execute(' INSERT INTO composition(CompID, Dict_Pop, pH) VALUES(?,?,?)', (CompID, compdictstr, pH))
    #         db.commit()
    #
    #
    #     except sqlite3.IntegrityError as e:
    #         print(str(e))
    #         if 'column CompID' in str(e):#.startswith('column OrgID'):
    #             print('\n CompID already in use, either update or use another \n')
    #             raise e
    #         elif 'UNIQUE' or 'columns' in str(e):#.startswith('columns'):
    #             print('\n This Composition already exists in the database \n')
    #
    #             cursor.execute('SELECT CompID FROM composition WHERE Dict_Pop = ? AND ' + \
    #                'pH = ?', (compdictstr, pH))
    #             ID = cursor.fetchone()
    #             print("\n Setting this composition's ID to " + ID[0] + "\n")
    #             self.CompID = ID[0]
    #             return ID[0]# set it here
    #         else:
    #             raise e
    #     finally:
    #         #ensure we safely close the database even if there is an error.
    #         db.close()

    def getconcs(self):
        """Returns the composition dicionary as a string."""
        compdictstr='{'
        for key in sorted(self.composition):
            compdictstr += "'" +key + "': " + str(self.composition[key].activity) +', '
            # compdict[key] = E.composition[key].activity
        compdictstr+='}'
        return compdictstr

    def update_pH(self, update, _from='pH'):
        """ Update the pH and [H+] of the reactor.

        Parameters
        ----------
        update : float
            The new value of pH or H+ concentration
        _from : str
            String identifier showing what ``update`` represents ('pH' or 'H+')
        """
        concH, pH = 0.,0.
        if _from=='pH':
            pH = update
            concH = 10**(-pH)
        elif _from=='H+':
            concH = update
            pH = - math.log10(concH)
        else:
            raise ValueError('Unclear what you are updating the pH with!')
        self.composition['H+'].conc = concH
        self.composition['H+'].activity = concH
        self.pH = pH



    # def update_locale_db(self, dbpath=nmp.std_dbpath):
    #     db = sqlite3.connect(str(dbpath))
    #     cursor = db.cursor()
    #
    #     try:
    #         cursor.execute(' INSERT INTO locale (LocID, EnvType) VALUES(?,?)',
    #           (self.LocID, self.tabname))
    #         db.commit()
    #     except sqlite3.IntegrityError as e:
    #         if 'column LocID' in str(e):
    #             print('\n LocID already in use, either update or use another \n')
    #             raise e
    #         else:
    #             raise e
    #     finally:
    #         db.close()
