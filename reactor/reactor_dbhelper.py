import sqlite3, os, ast
from datetime import date
import logging
from pandas import read_sql_query
from NutMEG.util.sqlshortcuts import sqlshortcuts as sqlsc
from NutMEG import reaction as rxn

import NutMEG.util.NutMEGparams as nmp
from NutMEG.util.loggersetup import loggersetup as logset
logger = logset.get_logger(__name__, filelevel=nmp.filelevel, printlevel=nmp.printlevel)

empty_default_dbdict = {'Pressure' : ['REAL', 0],
  'Temperature' : ['REAL', 0],
  'Volume' : ['REAL', 0],
  'pH' : ['REAL', 0],
  'reactions' : ['TEXT', ''],
  'CompID' : ['TEXT', ''],
  'composition_inputs' : ['TEXT', '{}']}

class rdb_helper:

    ReactIDs = []

    def __init__(self, host, dbpath=nmp.std_dbpath):

        self.host = host
        self.dbpath = dbpath

    def workoutID(self, dbdict=None):
        """Find this organisms' ID from the database, and if there's no
        entry create one. """

        db = sqlite3.connect(self.dbpath)
        cursor = db.cursor()

        self.update_ReactIDs()
        self.update_CompID()

        self.get_db_sqlparams(dbdict=dbdict)

        # search in the database for the organism's identifiers
        try:
            logger.info('Trying to find the reactor with name: ' + \
              self.host.name + ' in the database')
            command, vals = sqlsc.ANDlst(self.dbdict)
            cursor.execute(sqlsc.SELECTpreamble(self.host.name, 'LocID')+command, vals)
            op = cursor.fetchone()
            if op is not None:
                logger.info("Success! Its LocID is: "+op[0])
                self.host.LocID = op[0]
            else:
                raise Exception('No reactor found')
        except Exception as e:
            if 'no such table' in str(e):
                logger.info('Undefined table for this reactor, creating one.')
                self.createtable()
                logger.info("Table created. Adding reactor type: "+self.host.name+".")
                self.to_db()
            elif str(e) == 'No reactor found':
                logger.info('This is a unique reactor! Adding it to the ' + \
                  self.host.name + ' table.')
                self.to_db()
            else:
                logger.info("Other error encountered: "+str(e))
                logger.info('Assuming this is a new reactor of type: ' + \
                  self.host.name + '. Adding  to database.')
                logger.info('First make sure the table exists')
                self.createtable()
                logger.info('Now, adding an entry.')
                self.to_db()
        finally:
            db.close()

    ##### THIS REACTOR TYPE #########

    def createtable(self, replace=False):
        """Create a table in the database for this species."""
        db = sqlite3.connect(self.dbpath)
        cursor = db.cursor()

        if replace:
            cursor.execute('DROP TABLE: '+self.host.name)
            # also remove entries in the reactor table
            cursor.execute('DELETE FROM Reactor WHERE Type = '+self.host.name)


        ID_and_params = {'LocID' : ['TEXT', None]}
        ID_and_params.update(self.dbdict)

        cursor.execute('CREATE TABLE IF NOT EXISTS ' + self.host.name + ' (' +\
          sqlsc.commas_string_keys_types(ID_and_params) + ', PRIMARY KEY(LocID)' + \
          ', UNIQUE(' + sqlsc.commas_string_keys(self.dbdict) + '))')

        db.commit()
        db.close()


    def get_db_sqlparams(self, dbdict=None):
        """Get the parameters to import into the database as a dictionary"""
        if dbdict == None:
            self.dbdict = {'Pressure' : ['REAL', self.host.env.P],
              'Temperature' : ['REAL', self.host.env.T],
              'Volume' : ['REAL', self.host.env.V],
              'pH' : ['REAL', self.host.pH],
              'reactions' : ['TEXT', tuple(self.host.ReactIDs)],
              'CompID' : ['TEXT', self.host.CompID],
              'composition_inputs' : ['TEXT', str(self.host.composition_inputs)]}

        else:
            self.dbdict=dbdict


    def to_db(self):
        """Send the reactor's parameters to be stored in the database's table
        for this reactor type.

        If this eactor is already in the database, reset its LocID to match.
        """
        db = sqlite3.connect(self.dbpath)
        cursor = db.cursor()
        try:
            # attempt to add in an entry, with LocID 'Tester'
            Testerdict = {'LocID' : ['TEXT', 'Tester']}
            Testerdict.update(self.dbdict)

            command, vals = sqlsc.INSERTlst(self.host.name, Testerdict)
            cursor.execute(command, vals)
            db.commit()

            # if this goes through, generate an OrgID
            cursor.execute('SELECT rowid FROM ' + self.host.name + ' WHERE' + \
              ' LocID = ?', ('Tester',))
            entryno = cursor.fetchone()[0]
            self.host.LocID = (self.host.name + \
              str(entryno) + date.today().strftime("_%d%m%y"))


            # Now replace the 'Tester OrgID in the database'
            cursor.execute('UPDATE '+self.host.name + ' SET LocID = ? ' + \
              'WHERE LocID = ?', (self.host.LocID, 'Tester'))
            db.commit()

        except sqlite3.IntegrityError as e:
            if 'column LocID' in str(e):
                logger.INFO("\n LocID already in use, Tester hsan't been replaced somewhere \n")
                raise e
            elif 'columns' or 'UNIQUE' in str(e):
                logger.INFO('\n It looks like a reactor with these parameters ' + \
                  'already exists in the database \n')

                reset, vals = sqlsc.ANDlst(self.dbdict)
                cursor.execute(sqlsc.SELECTpreamble(self.host.name, 'LocID')+reset, vals)
                ID = cursor.fetchone()
                logger.INFO("\n Setting this Reactor's ID to " + ID[0] + "\n")
                self.host.LocID = ID[0]

            else:
                raise e
        finally:
            db.close()

        # Now the organism database needs to be updated so we can find these
        # details in ecosystem
        self.update_reactor_db()


    @staticmethod
    def from_db(name, LocID, dbdict=empty_default_dbdict,
      dbpath=nmp.std_dbpath):
        """ Extract the data from the database and return a dictionary
        of parameters to initialise the reactor.
        """

        db = sqlite3.connect(dbpath)
        cursor = db.cursor()
        # extract all the data in alphabetical order by column names
        cursor.execute(sqlsc.SELECTcolumns(name, dbdict, 'LocID'), (LocID,))
        Params = cursor.fetchone()
        db.close()

        for k, nv in zip(sorted(dbdict.keys()), Params):
            dbdict[k] = [dbdict[k][0], nv]

        return dbdict

    ####### reactor db ###########

    def update_reactor_db(self):
        """Update the Reactor table, which helps point to the reactor
        type tables."""
        db = sqlite3.connect(self.dbpath)
        cursor = db.cursor()
        #first, send the organism to the organism database
        try:
            cursor.execute(' INSERT INTO Reactor (LocID, Type) VALUES(?,?)',
              (self.host.LocID, self.host.name))
            db.commit()
        except sqlite3.IntegrityError as e:
            if 'column LocID' in str(e):
                logger.warning('LocID already in use, either update or ' + \
                ' use another')
                raise e
            else:
                raise e
        finally:
            db.close()

    ######## REACTIONS ######

    def update_ReactIDs(self):
        """If a new reaction has been added, use this to make sure all the
        ReactIDs match."""
        ReactIDlst = []
        for rxxneq in self.host.reactionlist.keys():
            for rxxn in self.host.reactionlist[rxxneq].values():
                ReactIDlst.append(self.add_to_Reactions(rxxn))
        self.host.ReactIDs = ReactIDlst

    @staticmethod
    def reaction_db_sqlparams(rxxn):
        """Get a dictionary of the values and datatypes to enter
        into the Reactions table. """
        reac, prod = rxxn.reagents_name()
        return {'Type' : ['TEXT', str(type(rxxn))],
          'equation' : ['TEXT', rxxn.equation],
          'reactants' : ['TEXT', str(reac)],
          'products' : ['TEXT', str(prod)]}

    def add_to_Reactions(self, rxxn):
        """Add the reaction passed to the Reactions table if it is not
        already saved. Return the ReactID of the reaction in the table."""

        db = sqlite3.connect(self.dbpath)
        cursor = db.cursor()
        try:
            # attempt to add in an entry, with ReactID 'Tester'
            Testerdict = {'ReactID' : ['TEXT', 'Tester']}
            rxxndict = rdb_helper.reaction_db_sqlparams(rxxn)
            Testerdict.update(rxxndict)

            command, vals = sqlsc.INSERTlst('Reactions', Testerdict)
            cursor.execute(command, vals)
            db.commit()

            # if this goes through, generate a ReactID
            cursor.execute('SELECT rowid FROM Reactions WHERE' + \
              ' ReactID = ?', ('Tester',))
            entryno = cursor.fetchone()[0]
            ReactID = str(entryno) + date.today().strftime("_%d%m%y")


            # Now replace the 'Tester OrgID in the database'
            cursor.execute('UPDATE Reactions SET ReactID = ? ' + \
              'WHERE ReactID = ?', (ReactID, 'Tester'))
            db.commit()

        except sqlite3.IntegrityError as e:
            if 'column ReactID' in str(e):
                logger.warning("ReactID already in use, Tester hsan't been replaced somewhere")
                raise e
            elif 'columns' or 'UNIQUE' in str(e):
                logger.info('It looks like a reaction with these ' + \
                  'parameters already exists in the database \n')

                reset, vals = sqlsc.ANDlst(rxxndict)
                cursor.execute(sqlsc.SELECTpreamble('Reactions', 'ReactID') + reset, vals)
                ID = cursor.fetchone()
                logger.info("Setting this Reaction's ID to " + ID[0])
                ReactID =  ID[0]
            else:
                raise e
        finally:
            db.close()

        return ReactID

    def extract_from_Reactions(self, ReactID, dbdict=None):
        """ Extract the data from the database and initialise a new reaction
        with its parameters.
        """

        # self.reaction_db_sqlparams(dbdict=dbdict)

        db = sqlite3.connect(self.dbpath)
        cursor = db.cursor()
        # extract all the data in alphabetical order by column names.
        # sloght fudge here because we only need the keys from
        # reaction_db_sqlparams
        tempdict = {'Type':['TEXT',0], 'equation':['TEXT',0],
          'reactants':['TEXT',0], 'products':['TEXT',0]}
        cursor.execute(sqlsc.SELECTcolumns('Reactions', tempdict, 'ReactID'), (ReactID,))
        Params = cursor.fetchone()
        db.close()

        for k, nv in zip(sorted(tempdict.keys()), Params):
            tempdict[k] = [tempdict[k][0], nv]

        # now rebuild a reaction from this data
        # concentrations are zero for now, once in reactor they will be
        # unified with its composition.
        reactants, products = {}, {}
        if tempdict['Type'][1] == str(rxn.reaction):
            for k, v in ast.literal_eval(tempdict['reactants'][1]).items():
                temp_reagent = rxn.reagent(k, self.host.env,
                  phase=rxn.reagent.get_phase_str(k))
                reactants[temp_reagent] = v

            for k, v in ast.literal_eval(tempdict['products'][1]).items():
                temp_reagent = rxn.reagent(k, self.host.env,
                  phase=rxn.reagent.get_phase_str(k))
                products[temp_reagent] = v

            return rxn.reaction(reactants, products, self.host.env)
        else:
            raise IntegrityError('NutMEG is not set up to recover this ' + \
              'type of reaction from the database yet!')


    ######### Composition ########

    # before we had comp_to_db and comp_from_db. Do something like reaction?

    def get_comp_activities(self):
        compdictstr='{'
        for key in sorted(self.host.composition):
            compdictstr += "'" + str(key) + "' : " + \
              str(self.host.composition[key].activity) +', '
        compdictstr+='}'
        return compdictstr


    def Composition_db_sqlparams(self):
        """Get a dictionary of the values and datatypes to enter
        into the Composition table. """
        return {'DictPop' : ['TEXT', self.get_comp_activities()],
          'pH' : ['REAL', self.host.pH]}


    def add_to_Composition(self):
        """Add the composition of the reactor to the Composition table if
        it is not already saved. Return the CompID."""

        db = sqlite3.connect(self.dbpath)
        cursor = db.cursor()

        compdictstr = self.get_comp_activities()

        try:
            # attempt to add in an entry, with CompID 'Tester'
            Testerdict = {'CompID' : ['TEXT', 'Tester']}
            compdict = self.Composition_db_sqlparams()
            Testerdict.update(compdict)

            command, vals = sqlsc.INSERTlst('Composition', Testerdict)
            cursor.execute(command, vals)
            db.commit()

            # if this goes through, generate a CompID
            cursor.execute('SELECT rowid FROM Composition WHERE' + \
              ' CompID = ?', ('Tester',))
            entryno = cursor.fetchone()[0]
            CompID = str(entryno) + date.today().strftime("_%d%m%y")

            # Now replace the 'Tester OrgID in the database'
            cursor.execute('UPDATE Composition SET CompID = ? ' + \
              'WHERE CompID = ?', (CompID, 'Tester'))
            db.commit()

        except sqlite3.IntegrityError as e:
            if 'column CompID' in str(e):
                logger.warning("CompID already in use, Tester hasn't ' + \
                  ' been replaced somewhere")
                raise e
            elif 'columns' or 'UNIQUE' in str(e):
                logger.info('It looks like a composition with these ' + \
                  'parameters already exists in the database \n')

                reset, vals = sqlsc.ANDlst(compdict)
                cursor.execute(sqlsc.SELECTpreamble('Composition', 'CompID') + reset, vals)
                ID = cursor.fetchone()
                logger.info("Setting this Composition's ID to " + ID[0])
                CompID =  ID[0]
            else:
                raise e
        finally:
            db.close()

        return CompID


    def extract_from_Composition(self, CompID):
        db = sqlite3.connect(self.dbpath)
        cursor = db.cursor()
        cursor.execute("SELECT * FROM Composition WHERE CompID = ?", (CompID,))

        params = cursor.fetchone()
        Comp = ast.literal_eval(params[1])
        pH = params[2]

        db.close()

        fullcomp = {}
        for k, v in Comp.items():
            temp_reagent = rxn.reagent(k, self.host.env,
              phase=rxn.reagent.get_phase_str(k), conc=v, activity=v)
            fullcomp[temp_reagent.name] = temp_reagent

        # include H+
        H = rxn.reagent('H+', self.host.env, charge=1, conc=10**-pH,
          phase='aq', activity=10**-pH)
        fullcomp[H.name] = H

        self.host.composition.update(fullcomp)


    def update_CompID(self):
        self.host.CompID = self.add_to_Composition()

    ########## Other handy functions




    def print_table(self):
        """Print out the table containing data for the current host species"""
        db=sqlite3.connect(self.dbpath)
        print(read_sql_query("SELECT * FROM " + self.host.name, db))
        db.close()
