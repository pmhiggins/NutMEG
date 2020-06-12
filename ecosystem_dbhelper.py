import NutMEG.ecosystem
import pandas as pd
import sqlite3
from datetime import date
from itertools import chain
import numpy as np
import os, ast, re
from NutMEG.util.sqlshortcuts import sqlshortcuts as sqlsc

import NutMEG.util.NutMEGparams as nmp
from NutMEG.util.loggersetup import loggersetup as logset
logger = logset.get_logger(__name__, filelevel=nmp.filelevel, printlevel=nmp.printlevel)

Sumdict = {'SimID':['TEXT',0], 'OrgIDs' : ['TEXT',0], 'LocID' : 'TEXT', 'OrgNums' : ['TEXT',0], 'PeakGR': ['TEXT',0], 'FinBM_vol' : ['TEXT',0], 'FinBM_cells_tot' : ['REAL',0], 'FinBM_vol_tot' : ['REAL',0], 'FinBM_kg_tot' : ['REAL',0]}

class db_helper:
    """Class for managing the NutMEG database and its contents regarding
    ecosystem and simulations.

    This class remains a mess while I decide whether we'll need the
    unique_params again or not. For now that is all commented out.

    Attributes
    ----------
    e : ecosystem
        The ecosystem object performing simulations to save.
    dbpath : str
        path to the SQL database. If not passed it will be set to the default
        (NutMEG_db in your working directory). The default can be changed in
        NutMEG/utils/NuMEGparams.
    OrgNames : tuple
        tuple of organism names to determine table names
    OrgIDs : tuple
        tuple of OrgIDs
    OrgNums : tuple
        tuple of the populations of each organism type.
    """
    uniqueOrgs = ["Methanogen"]
    uniqueReas = ["EnceladusVessel"]

    def __init__(self, e, dbpath):
        self.e = e
        self.getOrgdetails()
        self.dbpath = dbpath
        # self.ReaName = e.r.name

    def getOrgdetails(self):
        """Update initial details of the culture"""
        self.OrgNames=()
        self.OrgIDs=()
        self.OrgNums=()
        for org in self.e.c.all():
            self.OrgNames += (org.name,)
            self.OrgIDs += (org.OrgID,)
            self.OrgNums += (org.get_population(), )

    def buildresultsdict(self):
        """Return the basic entries to the results dictionary, it will be
        updated in the organsim and reactor db helpers.
        """
        # first add in the base dictionary
        d = {'Time':[], 'Composition': []}

        # now add in all organism-specific parameters in the database
        # for OrgName in self.OrgNames:
        #     if OrgName in self.uniqueOrgs:
        #         d = self.adduniquebits_dict(OrgName, d)
        return d

    def createtable(self):
        """ Create a table for this simulation.

        If the culture and reactor combination has already been performed,
        return False and the SimID of that simulation.
        If not, add a row to the Summary table for this simulation, and create a
        FullResults_Sim_ table for the full output.

        Please make sure that the OrgIDs and LocIDs are up-to-date before
        working out the SimID.
        """
        # First, check whether this sim already exists
        cs, ID = self.check_summary(notrace=False) # keep Tester in the db.
        if not cs:
            # if Summary has these orgs, no need to create - it already exists.
            return False, ID
        db=sqlite3.connect(self.dbpath)
        cursor = db.cursor()

        # Now, create a new table for this Simulation
        # to generate simID, we need to know the number of the next
        # entry in Summary. check_summary() made Tester for us in its place.
        try:
            cursor.execute('SELECT rowid FROM Summary WHERE SimID = ?', (ID,))
            entryno = cursor.fetchone()[0]

            self.SimID = str(entryno)+date.today().strftime("_%d%m%y")

            cursor.execute('UPDATE Summary SET SimID = ? WHERE SimID = ?',
              (self.SimID, ID))
            # now the summary table is ready to populate upon completion.

            # create the table for this simulation's results.
            cursor.execute('CREATE TABLE IF NOT EXISTS FullResults_Sim_' + \
              self.SimID + ' (Time REAL, Composition REAL)')
            # add in a column for each output of each colony/horde
            for simdata in self.e.c.full_output_datatypes():
                cursor.execute('ALTER TABLE FullResults_Sim_' + self.SimID + \
                  ' ADD COLUMN '+simdata[0]+' '+simdata[1]+';')
            db.commit()
        except:
            raise
        finally:
            db.close()

        # this was for unique bits in specefic species

        # for OrgName in self.OrgNames:
        #     if OrgName in self.uniqueOrgs:
        #         self.adduniquebits(OrgName)
        # if self.ReaName in self.uniqueReas:
        #     self.adduniquebits(self.ReaName)

        return True, self.SimID


    # ############### If we add back in unique bits  ###############
    #
    # def adduniquebits(self, bitname):#, curs, db):
    #     if bitname == 'Methanogen':
    #         # add in the methanogen bits
    #         self.addcolumn('FullResults_Sim_'+self.SimID, 'ThrottledGrowthPower TEXT')#, curs, db)
    #         self.addcolumn('FullResults_Sim_'+self.SimID, 'WastedCO2 TEXT')
    #
    # def adduniquebits_dict(self, bitname, dicti):#, curs, db):
    #     if bitname == 'Methanogen':
    #         # add in the methanogen bits
    #         dicti['ThrottledGrowthPower'] = np.array([])
    #         dicti['WastedCO2'] = np.array([])
    #     return dicti
    #
    # ###############  ###############

    def dict_to_db(self, dicti, end=False, finonly=False):
        """Save all info in dicti to the FullResults table for this simulation,
        and return the refreshed (emptied) dictionary.


        Parameters
        ----------
        dicti : dict
            dictionary to update FullResults and/or Summary with.
        end : bool, optional
            True if this is the final dict to be saved, then Summary is updated
            too. Default False
        finonly : bool, optional
            If True, only add in the last entry of dicti. Default False
        """

        # get the contents from culture_output
        dicti.update(self.e.c.full_output())
        db=sqlite3.connect(self.dbpath)
        cursor = db.cursor()


        names = []
        values = []
        for key in sorted(dicti):
            names.append(key)
            values.append(dicti[key])
        names = str(tuple(names))
        values = tuple(values)

        if not finonly:
            for i in range(0,len(dicti['Time'])):
                values = []
                for key in sorted(dicti):

                    values.append(dicti[key][i])
                qmarks='(?'
                for r in range(1,len(values)):
                    qmarks += ',?'
                qmarks += ')'
                values = tuple(values)

                cursor.execute(' INSERT INTO FullResults_Sim_'+self.SimID+
                  ' '+str(names)+' VALUES'+str(qmarks), values)
        else:
            values = []
            for key in sorted(dicti):
                values.append(dicti[key][-1])
            qmarks='(?'
            for r in range(1,len(values)):
                qmarks += ',?'
            qmarks += ')'
            values = tuple(values)

            cursor.execute(' INSERT INTO FullResults_Sim_'+self.SimID+
              ' '+str(names)+' VALUES'+str(qmarks), values)

        db.commit()
        db.close()

        if end:
            self.update_summary(dicti)

        self.e.c.refresh_output()

        return dicti


    def update_summary(self, dicti):
        """At the end of the simulation, update the Summary table using the
        final entry in dicti."""

        db=sqlite3.connect(self.dbpath)
        curs = db.cursor()
        try:

            curs.execute('UPDATE Summary SET FinBM_cells_tot = ? WHERE SimID = ?',
              (dicti['totBM_cells'][-1], self.SimID))
            curs.execute('UPDATE Summary SET FinBM_vol_tot = ? WHERE SimID = ?',
              (dicti['totBM_vol'][-1], self.SimID))
            curs.execute('UPDATE Summary SET FinBM_kg_tot = ? WHERE SimID = ?',
              (dicti['totBM_kg'][-1], self.SimID))

            PeakGR = tuple()
            vols = tuple()

            for n in self.OrgIDs:
                curs.execute('SELECT MAX(GrowthRate_'+n+') FROM' + \
                  ' FullResults_Sim_'+self.SimID)
                PeakGR += (curs.fetchone()[0],)
                vols += (dicti['Volume_'+n][-1],)

            curs.execute('UPDATE Summary SET PeakGR = ? WHERE' + \
              ' SimID = ?', (str(PeakGR), self.SimID))
            curs.execute('UPDATE Summary SET FinBM_vol = ? WHERE' + \
              ' SimID = ?', (str(vols), self.SimID))

            db.commit()
        except:
            raise
        finally:
            db.close()


    def addcolumn(self, Table, NameType):
        """Add a column to a table if it doesn't exist."""
        db=sqlite3.connect(self.dbpath)
        curs = db.cursor()
        try:
            curs.execute('ALTER TABLE '+Table+' ADD COLUMN '+NameType+';')
            db.commit()
        except:
            print('Looks like we already added the '+NameType+' column')
        finally:
            db.close()


    def check_summary(self,notrace=True):
        """Check to see if the simulation you have set up has been performed
        before.

        If it has, a False value and the ID of that simulaiton will be returned.
        If not, a True value and a Tester ID of a new entry in Summary will be
        returned.

        Parameters
        ----------
        notrace : bool, optional
            If True, delete the entry after checking. Default is True
        """
        db=sqlite3.connect(self.dbpath)
        curs = db.cursor()

        # make sure we have all the culture details, in case new
        # orgs/colines/hordes have been added
        self.getOrgdetails()

        try:

            curs.execute('INSERT INTO Summary(SimID, OrgIDs, LocID, OrgNums) VALUES(?,?,?,?)', ('Tester', str(self.OrgIDs), self.e.r.LocID, str(self.OrgNums)))
            if notrace:
                curs.execute('DELETE FROM Summary WHERE SimID = ?', ('Tester',))
            db.commit()
            return True, 'Tester'

        except sqlite3.IntegrityError as e:
            if 'columns' or 'UNIQUE' in str(e):
                logger.info('This combination of culture and location has already been performed')

                curs.execute('SELECT SimID from Summary WHERE OrgIDs = ? '+\
                  'AND LocID = ? AND OrgNums = ?', (str(self.OrgIDs), self.e.r.LocID, str(self.OrgNums)))
                self.SimID = curs.fetchone()[0]

                logger.info('  The SimID of this is: '+self.SimID)

                return False, self.SimID
            else:
                raise e
        finally:
            db.close()


    # def nogrowth(self):
    #     """Save es conditions if we can't find a timestep"""
    #
    #     resultsdict = self.buildresultsdict()
    #     startpop=[1]
    #     endpop=[1]
    #
    #     resultsdict['Time'] = np.append(resultsdict['Time'], str(0))
    #     resultsdict['no_alive'] = np.append(resultsdict['no_alive'], tuple(endpop))
    #     resultsdict['EnergyAvailable'] = np.append(resultsdict['EnergyAvailable'], tuple([-col.respiration.net_pathway.molar_gibbs for col in self.e.c.collection]))
    #     resultsdict['EnergyConservable'] = np.append(resultsdict['EnergyConservable'], tuple([col.respiration.G_C for col in self.e.c.collection]))
    #     #CO2molality.append(self.e.r.composition['CO2(aq)'].activity)
    #     resultsdict['CatabolicRatePerCell'] = np.append(resultsdict['CatabolicRatePerCell'], tuple([col.respiration.rate/sp for col, sp in zip(self.e.c.collection, startpop)]))
    #     resultsdict['MaintenanceFrac'] = np.append(resultsdict['MaintenanceFrac'], tuple(self.e.c.get_average_maintenance_fraction()))
    #     #no_cells.append(self.e.c.get_population())
    #     #consumedC.append(self.e.c.collection[0].respiration.rate*dt)
    #     resultsdict['ColonyVol'] = np.append(resultsdict['ColonyVol'], tuple([col.totvolume for col in self.e.c.collection]))
    #     resultsdict['totColonyVol'] = np.append(resultsdict['totColonyVol'], tuple(self.e.c.volume))
    #     ####
    #     resultsdict['PowerSupplyPerCell'] = np.append(resultsdict['PowerSupplyPerCell'], tuple([col.P_s/sp for col, sp in zip(self.e.c.collection, startpop)]))
    #     resultsdict['GrowthPowerPerCell'] =  np.append(resultsdict['GrowthPowerPerCell'], tuple([col.P_growth/sp for col, sp in zip(self.e.c.collection, startpop)]))
    #
    #     #previously dataheavy:
    #     resultsdict['GrowthRate'] = np.append(resultsdict['GrowthRate'], tuple([0 for ep, sp in zip(endpop, startpop)]))
    #     resultsdict['Composition'] = np.append(resultsdict['Composition'], self.e.r.getconcs())
    #     resultsdict['CHNOPSUptakes'] = np.append(resultsdict['CHNOPSUptakes'], tuple(self.e.c.getCHNOPSut()))
    #     resultsdict['Limiter'] = np.append(resultsdict['Limiter'], tuple([col.CHNOPS.limiter for col in self.e.c.collection]))
    #     resultsdict['Throttling'] = np.append(resultsdict['Throttling'], tuple([col.throttling for col in self.e.c.collection]))
    #     # timeenddict=timer.time()
    #     # print("Saving to dict time: "+ str(timeenddict-timestartdict))
    #
    #     # ake unique organism things the new dataheavy
    #     # if dataheavy:
    #     #     resultsdict = self.e.CRuniquemonitors(resultsdict, dt, startpop)
    #
    #     self.dict_to_db(resultsdict, end=True)


    @staticmethod
    def findOrgIDsLocID(SimID, dbpath=nmp.std_dbpath):
        """Return the OrgIDs and LocID from the simualation that had SimID"""
        db=sqlite3.connect(dbpath)
        curs = db.cursor()

        try:
            curs.execute('SELECT OrgIDs FROM Summary WHERE SimID = ?', (SimID,))
            OrgIDs = curs.fetchone()[0]
            curs.execute('SELECT LocID FROM Summary WHERE SimID = ?', (SimID,))
            LocID = curs.fetchone()[0]
            return OrgIDs, LocID
        except sqlite3.Error as e:
            print(str(e))
            print('\n This SimID could not be found, are you sure it has been performed? \n')

            return None
        finally:
            db.close()


    @staticmethod
    def findSimID(OrgIDs, LocID, OrgNums, dbpath=nmp.std_dbpath):
        """Look in the Summary table to see if a simulation has been done
        with these organisms and environments. If it has, return the SimID.

        Parameters
        ----------
        OrgIDs : tuple
            OrgIDs included in the simulation, in the correct order
        LocID : str
            LocID to extract.
        OrgNums : tuple
            Initial number of organism, in the correct order.
        dbpath : str, optional
            path to database. Default is NutMEG_db
        """

        if isinstance(OrgIDs, str):
            # if a string is passed there is only one organism.
            OrgIDs = (OrgIDs,)
            OrgNums = (OrgNums,)

        db=sqlite3.connect(dbpath)
        curs = db.cursor()

        try:
            curs.execute('SELECT SimID FROM Summary WHERE OrgIDs = ? AND LocID = ? AND OrgNums = ?', (str(OrgIDs), LocID, str(OrgNums)))
            SimID = curs.fetchone()[0]
            return SimID
        except sqlite3.Error as e:
            print(str(e))
            logger.warning('This combination of Organism(s) and location could not be found, are you sure it has been performed?')

            return None
        finally:
            db.close()



    @staticmethod
    def extract_param_db_Sim(SimID, param, OrgID='', dbpath=nmp.std_dbpath):
        """Fetch a specific parameter from a simulation, using its ID as
        passed. If the param is a Summary variable this returns the final
        value from that simulation, if its one of the others, return
        the full list of results.

        Parameters
        ----------
        SimID : str
            Simulation ID to use as identifier in Summay or FullResults tables.
        param : str
            Parameter name to extract. Should match database column name,
            equivalent to parameters found in horde_output, culture_output etc.
        OrgID : str, optional
            Some parameters are organism specific and so require the OrgID to
            extract. If this is the case, pass the OrgID. Default ''
        dbpath : str, optional
            path to database. Default is NutMEG_db
        """
        db=sqlite3.connect(dbpath)
        cursor = db.cursor()
        colnames = Sumdict.keys()

        try:
            for c in colnames:
                if param == c:
                    # this is in Summary, which has one values for the simulation.
                    cursor.execute('SELECT '+param+' FROM Summary WHERE SimID = ?', (SimID,))
                    return cursor.fetchone()[0]

            # if the param is not in Summary, look for it in FullResults.
            if OrgID=='':
                com= 'SELECT '+param+' FROM FullResults_Sim_'+SimID
            else:
                com= 'SELECT '+param+'_'+OrgID+' FROM FullResults_Sim_'+SimID
            cursor.execute(com)

            return cursor.fetchall()

        except sqlite3.Error or TypeError as e:
            print(e)
            print('Problem retrieving data, check SimID and requested variables')
            return [0]
        except Exception as e:
            raise e
        finally:
            db.close()


    @staticmethod
    def extract_param_db_OrgLoc(OrgIDs, LocID, OrgNums, param, dbpath=nmp.std_dbpath):
        """Fetch a specific parameter from a simulation, using the relevant
        OrgID and LocIDs as passed. If the param is a Summary variable, this
        returns the final value from that simulation., if its one of the
        others, this returns the full list of results

        TODO: add in support for organism-specific parameters

        Parameters
        ----------
        OrgIDs : tuple
            OrgIDs included in the simulation, in the correct order
        LocID : str
            LocID to extract.
        OrgNums : tuple
            Initial number of organism, in the correct order.
        param : str
            Parameter name to extract. Should match database column name,
            equivalent to parameters found in horde_output, culture_output etc.
        dbpath : str, optional
            path to database. Default is NutMEG_db
        """

        SimID = db_helper.findSimID(OrgIDs, LocID, OrgNums, dbpath=dbpath)
        return db_helper.extract_param_db_Sim(SimID, param, dbpath=dbpath)



    @staticmethod
    def removesimdata(SimID, dbpath=nmp.std_dbpath):
        """ Delete data from a simulation"""
        # TODO: add in an option to delete the org and loc data too
        db=sqlite3.connect(dbpath)
        cursor = db.cursor()
        # delete the entry from the summary table and
        try:

            cursor.execute('DROP TABLE FullResults_Sim_'+SimID)
            cursor.execute("DELETE FROM Summary WHERE SimID = ?", (SimID,))
            db.commit()
        except Exception as e:
            raise e
        finally:
            db.close()

    @staticmethod
    def print_table(tabname, dbpath=nmp.std_dbpath):
        """print out table tabname using pandas"""
        db=sqlite3.connect(dbpath)
        try:
            print('\t '+tabname+' \n')
            print(pd.read_sql_query("SELECT * FROM "+tabname+' ', db))
        except:
            raise
        finally:
            db.close()



    @staticmethod
    def print_major_tables(dbpath=nmp.std_dbpath, FullResults=None):
        """Print out the major databases (e.g. all without the sims).

        If you also want to pring out a FullResults table, pass
        FullResults as the SimID."""

        db=sqlite3.connect(dbpath)
        try:
            print('\t Reactor \n')
            print(pd.read_sql_query("SELECT * FROM Reactor ", db))
            print('\n \t Composition \n')
            print(pd.read_sql_query("SELECT * FROM Composition ", db))
            print('\n \t Reactions \n')
            print(pd.read_sql_query("SELECT * FROM Reactions", db))
            # print('\n \t Enceladus \n')
            # print(pd.read_sql_query("SELECT * FROM Enceladus ", db))
            print('\n \t Organism \n')
            print(pd.read_sql_query("SELECT * FROM Organism ", db))
            # print('\n \t Methanogen \n')
            # print(pd.read_sql_query("SELECT * FROM Methanogen ", db))
            print('\n \t Summary \n')
            print(pd.read_sql_query("SELECT * FROM Summary ", db))

            if FullResults != None:
                print('\n \t Full Results '+FullResults+ '\n')
                print(pd.read_sql_query("SELECT * FROM FullResults_Sim_"+FullResults, db))
        except:
            raise
        finally:
            db.close()



    @staticmethod
    def create_major_tables(replace=False, dbpath=nmp.std_dbpath):
        """Create the major tables in the ecosystem database inside the
        file passsed. If replace is True, any tables of the same name will
        be deleted. """

        db = sqlite3.connect(dbpath)
        cursor = db.cursor()
        if replace:
            cursor.execute('DROP TABLE Reactor')
            cursor.execute('DROP TABLE Composition')
            cursor.execute('DROP TABLE Reactions')
            # cursor.execute('DROP TABLE Enceladus')
            # cursor.execute('DROP TABLE Venus')
            cursor.execute('DROP TABLE Organism')
            # cursor.execute('DROP TABLE Methanogen')
            # cursor.execute('DROP TABLE SulfateReducer')
            cursor.execute('DROP TABLE Summary')

        cursor.execute('CREATE TABLE IF NOT EXISTS Reactor ' + \
          '(LocID TEXT, Type TEXT, PRIMARY KEY(LocID))')

        cursor.execute('CREATE TABLE IF NOT EXISTS Composition ' + \
          '(CompID TEXT, DictPop TEXT, pH REAL, PRIMARY KEY(CompID), ' + \
          'UNIQUE(DictPop, pH))')

        cursor.execute('CREATE TABLE IF NOT EXISTS Organism ' + \
          '(OrgID TEXT, Type TEXT, PRIMARY KEY(OrgID))')

        cursor.execute('CREATE TABLE IF NOT EXISTS Reactions' + \
          '(ReactID TEXT, Type TEXT, equation TEXT, ' + \
          'reactants TEXT, products TEXT, PRIMARY KEY(ReactID), ' + \
          'UNIQUE(Type, equation))')


        cursor.execute('CREATE TABLE IF NOT EXISTS Summary('+ sqlsc.commas_string_keys_types(Sumdict) + ', PRIMARY KEY(SimID), UNIQUE(OrgIDs, LocID, OrgNums))')

        # '(SimID TEXT, OrgIDs TEXT, LocID TEXT, OrgNums TEXT, ' + \
        # 'PeakGR TEXT, FinBM_vol TEXT, FinBM_cells_tot REAL, ' + \
        # 'FinBM_vol_tot REAL, FinBM_kg_tot REAL, ' + \


          # cursor.execute('CREATE TABLE IF NOT EXISTS Methanogen ' + \
          #   '(OrgID TEXT, Reducer TEXT, Esynth REAL, dry_mass REAL, ' + \
          #   'MP REAL DEFAULT 0, Volume REAL, Tdef TEXT, pHdef TEXT, ' + \
          #   'memb_pot REAL, PermH REAL, pHint REAL, Wasteful TEXT, ' + \
          #   'num INTEGER, n_ATP REAL, k_RTP REAL, PRIMARY KEY(OrgID), ' + \
          #   'UNIQUE(Reducer, Esynth, dry_mass, MP, Volume, Tdef, ' + \
          #   'pHdef, memb_pot, PermH, pHint, Wasteful, num, n_ATP, k_RTP))')
          #
          # cursor.execute('CREATE TABLE IF NOT EXISTS SulfateReducer ' + \
          #   '(OrgID TEXT, Reducer TEXT, Esynth REAL, dry_mass REAL, ' + \
          #   'MP REAL DEFAULT 0, Volume REAL, Tdef TEXT, pHdef TEXT, ' + \
          #   'memb_pot REAL, PermH REAL, pHint REAL, Wasteful TEXT, ' + \
          #   'num INTEGER, n_ATP REAL, k_RTP REAL, PRIMARY KEY(OrgID), ' + \
          #   'UNIQUE(Reducer, Esynth, dry_mass, MP, Volume, Tdef, ' + \
          #   'pHdef, memb_pot, PermH, pHint, Wasteful, num, n_ATP, k_RTP))')

        db.commit()
        db.close()


    @staticmethod
    def guess_name_from_ID(ID):
        """Try to guess an object's name from its ID"""
        nl = ID.split('_')
        nll = re.findall(r"[^\W\d_]+|\d+", str(nl))
        return nll[0]
