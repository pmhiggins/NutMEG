"""

Class for thermodynamic calculations regarding reactions.
This class isn't used much any more, but I'll leave it in
for legacy reasons.

@author P. M. Higgins

"""

from reaktoro import SupcrtDatabase
import math
import os


class reaction_thermo:

    def __init__(self, host, db="supcrt07-organics"):
        self.host=host
        self.dbname=db


        self.T, self.P = None, None
        self.get_stdG_lnK() # calculates dG and lnK in current rtr.env.

    def get_db(self):
        if self.dbname == None:
            _db = SupcrtDatabase("supcrt07-organics")
        elif type(self.dbname) == type(' '):
            _db = SupcrtDatabase(self.dbname)
        else:
            _db = self.dbname
        return _db


    def get_stdG_lnK(self, T=None, P=None):

        T, P = self.TPcheck(T,P)

        if self.T == T and self.P == P:
            # T and P have not changed, no need to recalculate
            return self.stdG, self.lnK
        else:
            # they have changed, recalculate.
            self.T = T
            self.P = P
            rkt_rxn = self.get_db().reaction(self.host.equation)
            rprops = rkt_rxn.props(self.T, 'K', self.P, 'Pa')
            self.lnK = rprops.lgK * math.log(10.)
            self.stdG = rprops.dG0
            return self.stdG, self.lnK

    def TPcheck(self, T, P):
        if T==None and P==None:
            return round(self.host.env.T,2), round(self.host.env.P,2)
        else:
            return round(T,2), round(P,2)
