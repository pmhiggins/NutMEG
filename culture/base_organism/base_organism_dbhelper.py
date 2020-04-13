# from NutMEG
import sqlite3, os
from datetime import date
import logging
from pandas import read_sql_query
from NutMEG.util.sqlshortcuts import sqlshortcuts as sqlsc

import NutMEG.util.NutMEGparams as nmp
from NutMEG.util.loggersetup import loggersetup as logset
logger = logset.get_logger(__name__, filelevel=nmp.filelevel, printlevel=nmp.printlevel)

empty_default_dbdict = {'Respiration' : ['TEXT', ''],
  'Esynth' : ['REAL', 0],
  'DryMass' : ['REAL', 0],
  'Mass' : ['REAL', 0],
  'MaintenancePower' : ['TEXT', ''],
  'Volume' : ['REAL', 0],
  'Tdef' : ['TEXT', ''],
  'pHdef' : ['TEXT', ''],
  'MembranePot' : ['REAL', 0],
  'PermH' : ['REAL', 0],
  'PermOH' : ['REAL', 0],
  'pHint' : ['REAL', 0],
  'n_ATP' : ['REAL', 0],
  'k_RTP' : ['REAL', 0],
  'base_life_span' : ['REAL', float('inf')]
  }

class bodb_helper:

    def __init__(self, host, dbpath=nmp.std_dbpath):

        self.host = host
        self.dbpath = dbpath

    def workoutID(self, dbdict=None):
        """Find this organisms' ID from the database, and if there's no
        entry create one. """

        db = sqlite3.connect(self.dbpath)
        cursor = db.cursor()

        self.get_db_sqlparams(dbdict=None)

        # search in the database for the organism's identifiers
        try:
            logger.debug('Trying to find your organism with name: ' + \
              self.host.name + ' in the database')
            command, vals = sqlsc.ANDlst(self.dbdict)
            cursor.execute(sqlsc.SELECTpreamble(self.host.name, 'OrgID')+command, vals)
            op = cursor.fetchone()
            if op is not None:
                logger.debug("Success! Its OrgID is: "+op[0])
                self.host.OrgID = op[0]
            else:
                raise Exception('No organism found')
        except Exception as e:
            if 'no such table' in str(e):
                logger.info('Undefined table for this organism, creating one.')
                self.createtable()
                logger.info("Table created. Adding species: "+self.host.name+".")
                self.to_db()
            elif str(e) == 'No organism found':
                logger.info('This is a unique organism! Adding it to the ' + \
                  self.host.name + ' table.')
                self.to_db()
            else:
                logger.info("Other error encountered: "+str(e))
                logger.info('Assuming this is a new organism of species: ' + \
                  self.host.name + '. Adding  to database.')
                self.to_db()
        finally:
            db.close()


    def createtable(self, replace=False):
        """Create a table in the database for this species."""
        db = sqlite3.connect(self.dbpath)
        cursor = db.cursor()

        if replace:
            cursor.execute('DROP TABLE: '+self.host.name)
            # also remove entries in the organism table
            cursor.execute('DELETE FROM organism WHERE Type = '+self.host.name)


        ID_and_params = {'OrgID' : ['TEXT', None]}
        ID_and_params.update(self.dbdict)

        cursor.execute('CREATE TABLE IF NOT EXISTS ' + self.host.name + ' (' +\
          sqlsc.commas_string_keys_types(ID_and_params) + ', PRIMARY KEY(OrgID)' + \
          ', UNIQUE(' + sqlsc.commas_string_keys(self.dbdict) + '))')

        db.commit()
        db.close()



    def get_db_sqlparams(self, dbdict=None):
        """Get the parameters to import into the database as a dictionary"""
        if dbdict == None:
            self.dbdict = {'Respiration' : ['TEXT', self.host.respiration.net_pathway.equation],
              'Esynth' : ['REAL', self.host.E_synth],
              'DryMass' : ['REAL', self.host.dry_mass],
              'Mass' : ['REAL', self.host.mass],
              'MaintenancePower' : ['TEXT', self.host.maintenance.get_netdictstr()],
              'Volume' : ['REAL', self.host.base_volume],
              'Tdef' : ['TEXT', self.host.maintenance.Tdef],
              'pHdef' : ['TEXT', self.host.maintenance.pHdef],
              'MembranePot' : ['REAL', self.host.memb_pot],
              'PermH' : ['REAL', self.host.PermH],
              'PermOH' : ['REAL', self.host.PermOH],
              'pHint' : ['REAL', self.host.pH_interior],
              'n_ATP' : ['REAL', self.host.respiration.n_ATP],
              'k_RTP' : ['REAL', self.host.respiration.net_pathway.rate_constant_RTP],
              'base_life_span' : ['REAL', self.host.base_life_span]
              }
        else:
            self.dbdict=dbdict

    @staticmethod
    def from_db(name, OrgID, dbdict=empty_default_dbdict, dbpath=nmp.std_dbpath):
        """ Extract the data from the database and return a dictionary
        of parameters.
        """

        db = sqlite3.connect(dbpath)
        cursor = db.cursor()
        # extract all the data in alphabetical order by column names
        cursor.execute(sqlsc.SELECTcolumns(name, dbdict, 'OrgID'), (OrgID,))
        Params = cursor.fetchone()
        db.close()

        for k, nv in zip(sorted(dbdict.keys()), Params):
            dbdict[k] = [dbdict[k][0], nv]

        return dbdict



    def to_db(self):
        """Send the organism's parameters to be stored in the database's table
        for this species.

        If this organism is already in the database, reset its OrgID to match.
        """

        db = sqlite3.connect(self.dbpath)
        cursor = db.cursor()
        try:
            # attempt to add in an entry, with OrgID 'Tester'
            Testerdict = {'OrgID' : ['TEXT', 'Tester']}
            Testerdict.update(self.dbdict)

            command, vals = sqlsc.INSERTlst(self.host.name, Testerdict)
            cursor.execute(command, vals)
            db.commit()

            # if this goes through, generate an OrgID
            cursor.execute('SELECT rowid FROM ' + self.host.name + ' WHERE' + \
              ' OrgID = ?', ('Tester',))
            entryno = cursor.fetchone()[0]
            self.host.OrgID = (self.host.name + \
              str(entryno) + date.today().strftime("_%d%m%y"))


            # Now replace the 'Tester OrgID in the database'
            cursor.execute('UPDATE '+self.host.name + ' SET OrgID = ? ' + \
              'WHERE OrgID = ?', (self.host.OrgID, 'Tester'))
            db.commit()

        except sqlite3.IntegrityError as e:
            if 'column OrgID' in str(e):
                logger.warning("\n OrgID already in use, Tester hsan't been replaced somewhere \n")
                raise e
            elif 'columns' or 'UNIQUE' in str(e):
                logger.info('\n It looks like an organism with these parameters ' + \
                  'already exists in the database \n')

                reset, vals = sqlsc.ANDlst(self.dbdict)
                cursor.execute(sqlsc.SELECTpreamble(self.host.name, 'OrgID')+reset, vals)
                ID = cursor.fetchone()
                logger.info("\n Setting this Organism's ID to " + ID[0] + "\n")
                self.host.OrgID = ID[0]

            else:
                raise e
        finally:
            db.close()

        # Now the organism database needs to be updated so we can find these
        # details in ecosystem
        self.update_org_db()


    def update_org_db(self):
        "Update the organism table, which points to the species tables."
        db = sqlite3.connect(self.dbpath)
        cursor = db.cursor()
        #first, send the organism to the organism database
        try:
            cursor.execute(' INSERT INTO Organism (OrgID, Type) VALUES(?,?)',
              (self.host.OrgID, self.host.name))
            db.commit()
        except sqlite3.IntegrityError as e:
            if 'column OrgID' in str(e):
                logger.warning('\n OrgID already in use, either update or use another \n')
                raise e
            else:
                raise e
        finally:
            db.close()


    def print_table(self):
        """Print out the table containing data for the current host species"""
        db=sqlite3.connect(self.dbpath)
        print(read_sql_query("SELECT * FROM " + self.host.name, db))
        db.close()




    ####### HELPFUL STRING GENERATING FUNCTIONS FOR SQLITE QUERIES ########
    #
    # def SELECTpreamble(self, SELECT='OrgID'):
    #     """Return the preamble to a SELECT Query"""
    #     return 'SELECT ' + SELECT + ' FROM ' + self.host.name + ' WHERE '
    #
    # def SELECTcolumns(self, dbdict, WHERE='OrgID'):
    #     """Return string for query to get all the values corresponding
    #     to the keys in dbdict"""
    #     First = True
    #     lst = 'SELECT '
    #     for k in sorted(dbdict.keys()):
    #         if First:
    #             First = False
    #             lst += k
    #             continue
    #         else:
    #             lst += ', '
    #             lst += k
    #     lst += ' FROM ' + self.host.name + ' WHERE ' + WHERE + ' = ?'
    #     return lst
    #
    #
    # def INSERTlst(self, dbdict):
    #     """Return string for query to insert all keys and values of dbdict"""
    #     First=True
    #     vals=[]
    #     lst = 'INSERT INTO ' + self.host.name + '(' + \
    #       self.commas_string_keys(dbdict) + ')'
    #     lst += ' VALUES('
    #     First=True
    #     for k in sorted(dbdict.keys()):
    #         vals.append(dbdict[k][1])
    #         if First:
    #             First = False
    #             lst += '?'
    #             continue
    #         else:
    #             lst += ','
    #             lst += '?'
    #     lst += ')'
    #
    #     return lst, tuple(vals)
    #
    # def commas_string_keys(self, dbdict):
    #     """Return keys of dbdict as one string separated by commas"""
    #     First=True
    #     lst= ''
    #     for k in sorted(dbdict.keys()):
    #         if First:
    #             First = False
    #             lst += k
    #             continue
    #         else:
    #             lst += ', '
    #             lst += k
    #     return lst
    #
    #
    # def commas_string_keys_types(self, dbdict):
    #     """Return keys of dbdict and their SQL data type as one string
    #     separated by commas"""
    #     First=True
    #     lst= ''
    #     for k in sorted(dbdict.keys()):
    #         if First:
    #             First = False
    #             lst += k + ' ' + dbdict[k][0]
    #             continue
    #         else:
    #             lst += ', '
    #             lst += k + ' ' + dbdict[k][0]
    #     return lst
    #
    #
    # def ANDlst(self, dbdict):
    #     """Return keys of dbdict in the sql AND x=? format. Also returns the
    #     corresponding vales as a tuple."""
    #     First = True
    #     vals=[]
    #     lst =''
    #     for k in sorted(dbdict):
    #         vals.append(dbdict[k][1])
    #         if First:
    #             First = False
    #             lst += k + ' = ?'
    #             continue
    #         else:
    #             lst += ' AND '
    #             lst += k + ' = ?'
    #     return lst, tuple(vals)
