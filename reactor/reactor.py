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
    """Module to store a reactive environment which we can name, such as
    methanogenesis at Enceladus.

    Holds an environment, and series of reactions in which we have an interest.
    Environmental parameters are in SI units unless otherwise specified.
    """

    env = environment()
    reactionlist = {} # dictionary of reactions of interest
    composition = {} #dictionary of reagents and their concentrations etc.
    volume = None
    pH = 7
    ReactIDs = tuple()

    def __init__(self, name, env=None, reactionlist={}, composition={},
      pH=7., workoutID=True, *args, **kwargs):
        self.name = name
        if env==None:
            self.env=environment()
        else:
            self.env = env
        self.reactionlist = reactionlist
        self.composition = composition
        self.volume = self.env.V
        self.CompID = ''
        self.pH=pH
        self.composition_inputs = kwargs.pop('composition_inputs', {})

        self.dbh= rdb_helper(self, dbpath=kwargs.pop('dbpath', nmp.std_dbpath))
        if workoutID:
            self.dbh.workoutID()

    @classmethod
    def r_from_db(cls, name, LocID, dbpath=nmp.std_dbpath):

        dbdict = rdb_helper.from_db(name, LocID, dbpath=dbpath)
        R = cls(name, env=environment(T = dbdict['Temperature'][1],
          P=dbdict['Pressure'][1], V=dbdict['Volume'][1]), pH=dbdict['pH'][1],
          workoutID=False, composition_inputs=ast.literal_eval(dbdict['composition_inputs'][1]))

        R.dbh.extract_from_Composition(dbdict['CompID'][1])
        R.rlist_from_ReactIDs(ast.literal_eval(dbdict['reactions'][1]))
        R.CompID = dbdict['CompID'][1]
        R.LocID = LocID
        return R


    def rlist_from_ReactIDs(self, ReactIDs):
        for rID in ReactIDs:
            self.add_reaction(self.dbh.extract_from_Reactions(rID))


    def print_reactions(self):
        """ print the list of equations, each one on a new line.
        """
        print(self.reactionlist.keys())

    def print_composition(self):
        print(self.composition.keys())

    def add_reaction(self, rxxn, overwrite=False):
        """Add a new reaction to reactionlist, unifying it with the
        environment.

        Pass overwrite as True if you want to overwrite the data we
        currently have about the reagents with what you have passed.
        """
        # add the reaction into reationlist
        self.reactionlist[rxxn.equation] = {type(rxxn):rxxn}
        self.unify_reaction(rxxn, overwrite=overwrite)

    def update_composition(self, t):
        """ If there are inflows into the composition, make the changes there
        would be in time t"""
        for c, r in self.composition_inputs.items():
            self.composition[c].activity += r*t
            self.composition[c].conc += r*t



    def unify_reaction(self, rxxn, overwrite=False):
        """
        Search the reactor for the reagents you've passed.
        If they exist here, unify them.
        If they don't, add them.
        The reaction passed will end up pointing to the composition of the
        reactor. Pass overwrite as True to update the values in the composition,
        False to use the values we already have.
        """
        #rrxnno = 0 # counting through the chain
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

        if overwrite:
            # overwrite the entry in reactionlist with the passed reacrion,
            # compositions and all.
            self.reactionlist[rxxn.equation][type(rxxn)] = rxxn

        # the reaction is not in the reactoinlist for this reactor, we
        # need to add it. First though, we must redefine the reagents
        # according to the composition.
        for rrxn in rxxn.reactants.keys():
            inlist = False
            for c_name, c_rxt in self.composition.items():
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



        for rrxn in rxxn.products.keys():
            inlist = False
            for c_name, c_rxt in self.composition.items():
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



    def perform_reaction(self, re_eq, n, re_type=rxn.reaction):
        """Perform the reaction with equation re_eq in reactionlist.

        n is the total number of unit molar reactions we want to take place.
        """
        self.reactionlist[re_eq][re_type].react(n)
        # update the composition at the end of the timestep in colony(lite)


    def change_vol(self, V):
        self.volume = float(V)
        self.env.V = float(V)

    def change_T(self, T):
        self.env.T = float(T)

    def change_P(self, P):
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
        compdictstr='{'
        for key in sorted(self.composition):
            compdictstr += "'" +key + "': " + str(self.composition[key].activity) +', '
            # compdict[key] = E.composition[key].activity
        compdictstr+='}'
        return compdictstr

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
