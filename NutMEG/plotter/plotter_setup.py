import matplotlib.pyplot as plt
import sqlite3
import os
import NutMEG
import NutMEG.util.NutMEGparams as nmp

std_dbpath=os.path.join(os.path.dirname(__file__), '../testDB')

class nutfig:
    """ Class for setting up plotting in NutMEG. Will end up doing
    data retrieval as well using the ecosystem_dbhelper"""

    fig = None

    def __init__(self, fig=None, figsize=(5,7)):
        if fig==None:
            self.fig = plt.figure(figsize=figsize)
        else:
            self.fig = fig

    @staticmethod
    def printcheck():
        """Checker method to see summary table. Probably isn't needed because
        you could just use ecosystem_dbhelper"""

        db=sqlite3.connect(std_dbpath)
        try:
            cursor = db.cursor()
            cursor.execute('SELECT * FROM Summary')
            rows = cursor.fetchall()

            for row in rows:
                print(row)
        except sqlite3.Error as e:
            raise e
        finally:
            db.close()


    @staticmethod
    def extract_param_db_Sim(SimID, param, **kwargs):
        """Fetch a specific parameter from a simulation, using its ID as
        passed. If the param is 'FinBM' or 'PeakGR' this returns the final
        value from that simulation, if its one of the others, this returns
        the full list of results"""
        return NutMEG.ecosystem_dbhelper.db_helper.extract_param_db_Sim(SimID, param, dbpath=kwargs.pop('dbpath', nmp.std_dbpath))
        # db=sqlite3.connect(std_dbpath)
        # db.row_factory = lambda cursor, row: row[0]
        # cursor = db.cursor()
        # try:
        #     if param == 'FinBM' or param == 'PeakGR':
        #         #use Summary, which has one value for the sims
        #         cursor.execute('SELECT '+param+' FROM Summary WHERE SimID = ?' (SimID))
        #         return cursor.fetchone()[0]
        #     else:
        #         com= 'SELECT '+param+' FROM FullResults_Sim_'+SimID
        #         cursor.execute(com)
        #
        #         return cursor.fetchall()
        # except:# sqlite3.Error or TypeError as e:
        #     print('Problem retrieving data, check SimID and requested variables')
        #     return [0]
        # finally:
        #     db.close()

    @staticmethod
    def extract_param_db_OrgLoc(OrgIDs, LocID, OrgNums, param, **kwargs):
        """Fetch a specific parameter from a simulation, using the relevant
        OrgID and LocIDs as passed. If the param is 'FinBM' or 'PeakGR' this returns the final value from that simulation OR if the simulation
        failed, a volume of 1e-18 and growth rate of 0 is returned. If its
        one of the others, this returns the full list of results."""
        return NutMEG.ecosystem_dbhelper.db_helper.extract_param_db_OrgLoc(OrgIDs, LocID, OrgNums, param, dbpath=kwargs.pop('dbpath', nmp.std_dbpath))
        # db=sqlite3.connect(std_dbpath)
        # db.row_factory = lambda cursor, row: row[0]
        # cursor = db.cursor()
        # try:
        #     if param == 'FinBM' or param == 'PeakGR':
        #         #use Summary, which has one value for the sims
        #         cursor.execute('SELECT '+param+' FROM Summary WHERE OrgIDs = ? AND LocID = ?', (str((OrgID, )), LocID))
        #         ans =cursor.fetchone()
        #         if ans==None and param=='FinBM':
        #             return 1e-18
        #         elif ans == None and param=='PeakGR':
        #             return 0
        #         else:
        #             return ans
        #     else:
        #         #Get the SimID to find the FullResults table
        #         cursor.execute('SELECT SimID FROM Summary WHERE OrgIDs = ? AND LocID = ?', (str((OrgID, )), LocID))
        #         SimID = cursor.fetchone()
        #         return nutfig.extract_param_db_Sim(SimID, param)
        # except sqlite3.Error as e:
        #     print('Problem retrieving data, make sure this Location and organism have been simulated together. Also check requested variables')
        #     return [0]
        #     # raise e
        # finally:
        #     db.close()

    @staticmethod
    def normalise(vals):
        """Take the list vals, and return it normalised between 0 and 1."""
        mi =min(vals)
        ma = max(vals)
        vr = []
        for v in vals:
            vr.append((v-mi)/(ma-mi))
        return vr, ma, mi

    @staticmethod
    def denormalise(val, maxval, minval):
        """Take a value val which sho be 'denormalised' from between
        0 and 1 to between maxval and minval."""
        return((val*(maxval-minval))+minval)
