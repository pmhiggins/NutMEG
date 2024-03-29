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
  'base_life_span' : ['REAL', float('inf')],
  'horde_deathrate' : ['REAL', 0.],
  'horde_num' : ['REAL', 500]
  }

class bodb_helper:
    """Class from managing how organisms and their attributes are stored in
    NutMEG's database.  Each organism/horde with a unique ``name``
    will have its own table, in case child classes of either have their
    own attributes they want saving. There is also an *organism* table, which stores
    all OrgIDs as pointers to the correct table.

    Attributes
    ----------
    host : base_organism like
        host organism/horde to be saving.
    dbpath : str
        path to the SQL database. If not passed it will be set to the default
        (NutMEG_db in your working directory). The default can be changed in
        NutMEG/utils/NuMEGparams.
    dbdict : dict
        Parameters which will be imported into the database.
    """
    def __init__(self, host, dbpath=nmp.std_dbpath):

        self.host = host
        self.dbpath = dbpath

    def get_db_sqlparams(self, dbdict=None):
        """Get the parameters to import into the database as a dictionary
        Pass a dbdict as a dictionary if you want to write your own."""
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
            try:
                self.dbdict['horde_deathrate'] = ['REAL', self.host.deathrate]
                self.dbdict['horde_num'] = ['REAL', self.host.num]
            except:
                # not a horde, should be fine though...
                self.dbdict['horde_deathrate'] = ['REAL', 0.]
                self.dbdict['horde_num'] = ['REAL', 1]
        else:
            self.dbdict=dbdict


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
            self.host.output.refreshparams()


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


    @staticmethod
    def extract_param_db(orgname, OrgID, params, dbpath=nmp.std_dbpath):
        """Fetch a specific parameter from a simulation, using its ID as
        passed."""
        db=sqlite3.connect(dbpath)
        cursor = db.cursor()
        try:
            sel = sqlsc.SELECTcolumns(orgname, params, 'OrgID')
            cursor.execute(sel, (OrgID,))
            return cursor.fetchall()

        except sqlite3.Error or TypeError as e:
            print(e)
            print('Problem retrieving data, check OrgID and requested variables')
            return [0]
        except Exception as e:
            raise e
        finally:
            db.close()
