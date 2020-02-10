"""

Class for thermodynamic calculations regarding reagents

"""

from reaktoro import Database, Thermo
import sqlite3 as sql
import os
import logging

db = Database("supcrt07-organics.xml")
thermo = Thermo(db)
rtodbpath = os.path.join(os.path.dirname(__file__), '../../data/TPdb')
# TPdb = sql.connect(dbpath)

import NutMEG.util.NutMEGparams as nmp
from NutMEG.util.loggersetup import loggersetup as logset
logger = logset.get_logger(__name__, filelevel=nmp.filelevel, printlevel=nmp.printlevel)


class reagent_thermo:

    def __init__(self, host):
        self.host=host


    def get_thermo_params(self, T=None, P=None):
        """ Use reaktoro to extract thermodynamic parameters for the host
        reagent.
        """
        T, P = self.TPcheck(T, P)
        species = self.host.name
        G = thermo.standardPartialMolarGibbsEnergy(T, P, species).val
        Hz= thermo.standardPartialMolarHelmholtzEnergy(T, P, species).val
        H = thermo.standardPartialMolarEnthalpy(T, P, species).val
        S = thermo.standardPartialMolarEntropy(T, P, species).val
        Cp= thermo.standardPartialMolarHeatCapacityConstP(T, P, species).val
        Cv= thermo.standardPartialMolarHeatCapacityConstV(T, P, species).val

        return G, Hz, H, S, Cp, Cv

    def get_RTP_params(self):
        """return the thermodynamic parameters at RTP. """
        return self.get_thermo_params(T=298.15, P=101325.0)

    def thermo_to_db(self, T=None, P=None):
        """Get ther thermodynamic parameters and save them to an SQLite
        database for retrieval elsewhere.
        """
        T, P = self.TPcheck(T, P)
        TPdb = sql.connect(rtodbpath)
        cursor = TPdb.cursor()
        try:
            cursor.execute('CREATE TABLE IF NOT EXISTS [' + self.host.name + \
              '] (T REAL, P REAL, G REAL, Hz REAL, ' + \
              'H REAL, S REAL, Cp REAL, Cv REAL)')

            G, Hz, H, S, Cp, Cv = self.get_thermo_params(T, P)
            cursor.execute(' INSERT INTO [' + self.host.name + ']' + \
              '(T, P, G, Hz, H, S, Cp, Cv) VALUES(?,?,?,?,?,?,?, ?)', (T, P,
              G, Hz, H, S, Cp, Cv))
            TPdb.commit()
            return True
        except:
            logger.warning('There is no thermodynamic data available for ' + \
            'your species: ' + self.host.name)
            return False
            #raise
        finally:
            TPdb.close()


    def db_select(self, paramschain='G, H, S, Cp', T=None, P=None):
        """Extract paramschain from the database for the host species"""
        T, P = self.TPcheck(T, P)
        TPdb = sql.connect(rtodbpath)
        cursor = TPdb.cursor()
        try:
            cursor.execute('SELECT ' + paramschain + ' FROM ['+self.host.name+
              '] WHERE T=? AND P=?', (T, P))
            data = cursor.fetchone()
            return data
        except:
            raise
        finally:
            TPdb.close()

    def TPcheck(self, T, P):
        """Ensure passed temperature and pressure is rounded to avoid
        crowding the database"""
        if T==None and P==None:
            return round(self.host.env.T,2), round(self.host.env.P,2)
        else:
            return round(T,2), round(P,2)
