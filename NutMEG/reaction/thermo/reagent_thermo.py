"""

Class for thermodynamic calculations regarding reagents

"""

from reaktoro import SupcrtDatabase
import os


import NutMEG.util.NutMEGparams as nmp
from NutMEG.util.loggersetup import loggersetup as logset
logger = logset.get_logger(__name__, filelevel=nmp.filelevel, printlevel=nmp.printlevel)


class reagent_thermo:

    def __init__(self, host, db="supcrt07-organics"):
        self.host=host
        self.dbname=db

        self.roundedT, self.roundedP = None, None
        #self.get_thermo_params() # in current rtr.env.


    def get_db(self):
        if self.dbname == None:
            _db = SupcrtDatabase("supcrt07-organics")
        elif type(self.dbname) == type(' '):
            _db = SupcrtDatabase(self.dbname)
        else:
            _db = self.dbname
        return _db


    def get_thermo_params(self, T=None, P=None, RTP=False):
        """ Use reaktoro to extract thermodynamic parameters for the host
        reagent.
        """
        _T, _P = self.TPcheck(T, P)

        if self.roundedT == _T and self.roundedP == _P:
            # T and P have not changed, no need to recalculate
            return self.G, self.Hz, self.H, self.S, self.Cp, self.Cv
        else:
            # they have changed or are not yet set, recalculate.

            # create a reaktoro species to extract the thermo information
            rkt_rxn = self.get_db().species(self.host.name)
            rprops = rkt_rxn.props(_T, 'K', _P, 'Pa')

            self.G = rprops.G0
            self.Hz= rprops.A0
            self.H = rprops.H0
            self.S = rprops.S0
            self.Cp= rprops.Cp0
            self.Cv= rprops.Cv0

            self.roundedT = _T
            self.roundedP = _P
            return self.G, self.Hz, self.H, self.S, self.Cp, self.Cv

    def get_RTP_params(self):
        """return the thermodynamic parameters at RTP. """
        return self.get_thermo_params(T=298.15, P=101325.0)


    def TPcheck(self, T, P):
        """Ensure passed temperature and pressure is rounded to avoid
        crowding the database"""
        if T==None and P==None:
            return round(self.host.env.T,2), round(self.host.env.P,2)
        else:
            return round(T,2), round(P,2)
