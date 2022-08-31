"""Module to store a reactive environment which we can name, such as
methanogenesis at Enceladus.

Holds an environment, composition and series of reactions in which
we have an interest.
Environmental parameters are in SI units unless otherwise specified.

Most recent changes: database fixes May 2020.

@author P M Higgins
"""
from NutMEG.environment import environment
from NutMEG import reaction as rxn
from itertools import chain
from copy import copy, deepcopy
import sqlite3
import sys, os, ast
from datetime import date
import time as timer
from .reactor_dbhelper import rdb_helper

import NutMEG.util.NutMEGparams as nmp
from NutMEG.util.loggersetup import loggersetup as logset
logger = logset.get_logger(__name__, filelevel=nmp.filelevel, printlevel=nmp.printlevel)

class reactor:
    """Class for storing reagents and reactions, and able to perform them.

    Attributes
    ----------
    name : str
        name of reactor. Used for database management (table name and LocID)
    reactionlist : dict
        reactions which may be performed in the format {'reaction eq' :
        {type(reaction) : NutMEG.reaction.reaction like}}. This allows for the
        same reaction to be included twice if we'd like to interpret it in
        different ways e.g. as a thermodynamic interaction or as a redox
        reaction.
    composition : dict
        reagents which can be found in the reactor in the format {
        'reagent name' : NutMEG.reaction.reagent like}
    pH : float, optional
        pH of the reactor, important for some interactions. Default is 7.0
    ReactIDs : tuple
        IDs of each reaction in the reactor, for saving to the database.
    env : environment, optional
        Local temperature, pressure and volume. Default RTP.
    composition_inputs : dict, kwarg
        If there is a net flow of reagents in/out of the reactor, set add them
        to this dictionary with their name as the key, and rate in M/s as the
        value.
    """

    def __init__(self, name, env=None, reactionlist={}, composition={},
      pH=7., workoutID=True, *args, **kwargs):
        self.name = name
        if env==None:
            self.env=environment()
        else:
            self.env = env
        self.reactionlist = reactionlist
        self.composition = composition
        self.volume = kwargs.pop('volume', self.env.V)
        self.env.V = self.volume
        self.CompID = ''
        self.pH=pH
        self.composition_inputs = kwargs.pop('composition_inputs', {})
        self.ReactIDs = tuple()

        self.dbh= rdb_helper(self, dbpath=kwargs.pop('dbpath', nmp.std_dbpath))
        if workoutID:
            self.dbh.workoutID()

    @classmethod
    def r_from_db(cls, name, LocID, dbpath=nmp.std_dbpath):
        """Extract a reactor from the SQL database at dbpath.

        Returns the saved reactor object

        Parameters
        ----------
        name : str
            name of the reactor. Required for table name.
        LocID : str
            LocID of the reaactor to extract.
        dbpath : str, optional
            location of the database file. Default is NutMEG_db outside the
            module directory.
        """

        dbdict = rdb_helper.from_db(name, LocID, dbpath=dbpath)
        R = cls(name, env=environment(T = dbdict['Temperature'][1],
          P=dbdict['Pressure'][1], V=dbdict['Volume'][1]), pH=dbdict['pH'][1],
          workoutID=False, composition_inputs=ast.literal_eval(dbdict['composition_inputs'][1]), dbpath=dbpath)

        R.dbh.extract_from_Composition(dbdict['CompID'][1])
        R.rlist_from_ReactIDs(ast.literal_eval(dbdict['reactions'][1]))
        R.CompID = dbdict['CompID'][1]
        R.LocID = LocID
        return R


    def rlist_from_ReactIDs(self, ReactIDs):
        """Set reactionlist from a list of ReactIDs"""
        for rID in ReactIDs:
            self.add_reaction(self.dbh.extract_from_Reactions(rID))


    def print_reactions(self):
        """ print the list of reaction equations.
        """
        print(self.reactionlist.keys())

    def print_composition(self):
        """print a list of the composition"""
        print(self.composition.keys())

    def add_reaction(self, rxxn, overwrite=False):
        """Add a new reaction rxxn to reactionlist, unifying it with the
        environment.

        Pass overwrite as True if you want to overwrite the data we
        currently have about the reagents with what you have passed.

        NB this function does not add new reagents to the reactor.
        """
        # add the reaction into reationlist
        self.reactionlist[rxxn.equation] = {type(rxxn):rxxn}
        self.unify_reaction(rxxn, overwrite=overwrite)

    def add_reagent(self, rct):
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



    def unify_reaction(self, rxxn, overwrite=False):
        """Add the reaction and its reagents to the reactor, ensuring there is
        only one of each reagent type in the reactor.

        Returns the reaction after unification

        Parameters
        ----------
        rxxn : NutMEG.reaction.reaction like
            reaction to unify
        overwrite : bool
            if True, overwrite the current composition with the activities of
            the reagents in rxxn. Default is False.

        Notes
        -----
        The reaction passed will end up pointing to the composition of the
        reactor. Pass overwrite as True to update the values in the composition,
        False to use the values we already have. If overwrite is False, and the
        reaction has a reagent which is not in this reactor's composition, it
        will be added.
        """

        if not overwrite:
            # we want the reaction passed to be reset to be the reaction
            # in reactionlist, if it exists in there.
            # if not, we need to add it, and unify the reagents with this
            # reactor's composition.
            try:
                # see if the reaction is already saved in the reactor
                rxxn = self.reactionlist[rxxn.name][type(rxxn)]
                # if it is, return that reaction
                return rxxn
            except:
                logger.info('Reaction to be unified not found in '+self.name)
                # return rxxn

        if overwrite:
            # overwrite the entry in reactionlist with the passed reacrion,
            # compositions and all.
            self.reactionlist[rxxn.equation][type(rxxn)] = rxxn

        # redefine the reagents
        # according to the composition.
        for rrxn in list(rxxn.reactants.keys()):
            inlist = False
            for c_name, c_rxt in self.composition.copy().items():
                inlist=True
                if rrxn.name == c_name:
                    # this reagent is in both the passed reaction and
                    # the composition.
                    if overwrite:
                        # update the composition entry with ragent data from
                        # the reaction
                        self.composition[c_name].redefine(rrxn)
                    # reset the reagent to be the same object as the one
                    # in the composition
                    rxxn.reactants[c_rxt] = rxxn.reactants.pop(rrxn)
            if not inlist:
                # it wasn't found in the composition, add it.
                logger.info('Adding '+rrxn.name+' to '+self.name+\
                  "'s composition.'")
                self.composition[rrxn.name] = rrxn

        for rrxn in list(rxxn.products.keys()):
            inlist = False
            for c_name, c_rxt in self.composition.copy().items():
                inlist=True
                if rrxn.name == c_name:
                    # this reagent is in both the passed reaction and
                    # the composition.
                    if overwrite:
                        # update the composition entry with ragent data from
                        # the reaction
                        self.composition[c_name].redefine(rrxn)
                    # reset the reagent to be the same object as the one
                    # in the composition
                    rxxn.products[c_rxt] = rxxn.products.pop(rrxn)

            if not inlist:
                # it wasn't found in the composition, add it.
                logger.info('Adding '+rrxn.name+' to '+self.name+\
                  "'s composition.'")
                self.composition[rrxn.name] = rrxn
        # return rxxn



    def perform_reaction(self, re_eq, n, re_type=rxn.reaction):
        """Perform the reaction n molar times.

        Parameters
        ----------
        re_eq : str
            reaction equation as it appears in reactionlist. If this equation is
            not in reactionlist an error is raised.
        n : float
            number of moles of reaction to perform.
        re_type : reaction like, optional
            The type of reaction to perform, if it is a subtype of
            NutMEG.reaction.reaction. Default is NutMEG.reaction.reaction.

        """
        self.reactionlist[re_eq][re_type].react(n)
        # update the composition at the end of the timestep in colony(lite)


    def change_vol(self, V):
        """Update reactor volume V in m^3"""
        self.volume = float(V)
        self.env.V = float(V)

    def change_T(self, T):
        """Update reactor temperature in K"""
        self.env.T = float(T)

    def change_P(self, P):
        """Update reactor pressure in Pa"""
        self.env.P = float(P)


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
        """ Update the pH of the reactor and [H+] to the float update. pass `from' as 'pH' to pass a pH value or 'H+' to pass a H+ concentration.
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
