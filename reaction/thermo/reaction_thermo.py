"""

Class for thermodynamic calculations regarding reactions.
This class isn't used much any more, but I'll leave it in
for legacy reasons.

@author P. M. Higgins

"""

from reaktoro import Database, Thermo
import sqlite3 as sql
import os

db = Database("supcrt07-organics.xml")
thermo = Thermo(db)
rtodbpath = os.path.join(os.path.dirname(__file__), '../../data/TPdb')


class reaction_thermo:

    def __init__(self, host):
        self.host=host


    def get_stdG_lnK(self, T=None, P=None):
        T, P = self.TPcheck(T,P)

        lnK = thermo.lnEquilibriumConstant(T, P, self.host.equation).val
        stdG = -8.314472*T*lnK

        return stdG, lnK

    def TPcheck(self, T, P):
        if T==None and P==None:
            return round(self.host.env.T,2), round(self.host.env.P,2)
        else:
            return round(T,2), round(P,2)
