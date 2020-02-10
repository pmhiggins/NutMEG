import sys
sys.path.append('..')
import NutMEG.reaction as reaction
from NutMEG.reaction.special.solutions import redox_half
from NutMEG.reaction.special.solutions import redox

import math
from NutMEG.environment import environment
from NutMEG.reactor import reactor
from itertools import chain
import ast
import sqlite3
import sys, os
from datetime import date

std_dbpath=os.path.join(os.path.dirname(__file__), '../../testDB')

R=0.08206 # L atm / mol K
gases = {'H2': {'SolubilityCons':[-39.9611,53.9381,16.3135,-0.036249,0.017565,-0.0023010], 'SolubilityForm':'Beta'},
  'CO2' : {'SolubilityCons':[-58.0931,90.5069,22.2940,0.027766,-0.025888,0.0050578], 'SolubilityForm':'K'}}

class VenusDrop(reactor):

    """The chemically active evirionment of the Enceladean Ocean.

    This is a saved instance of reative_system, which contains the
    methanogenesis and envrironmental data for the moon. In the full
    model this will be fleshed out with other constituents and maybe
    even other reactions/flows.
    """

    volume = 7.54e15 #m^3, from radius = 250km, depth = 10km
    env = environment(T=300., P=70e5, V=volume)
    depth=8.5
    LocID=''


    def __init__(self, H2ppm=25, HSact=1e-10, T=298, P=100000, setupcomp=True, autoID=True):
        self.env.T=298
        self.env.P=P
        self.name='VenusDrop'
        self.tabname = 'Venus'

        if setupcomp:
            self.initial_conditions(H2ppm, HSact,
              setupcomp=setupcomp)
            reactor.__init__(self, env=self.env,
              reactionlist=self.reactionlist,
              composition=self.composition)

        if autoID:
            self.workoutID()

    @classmethod
    def from_db(cls, LocID, dbpath=std_dbpath):
        """Load an Enceladus from the SQL database, table name: Enceladus
        Pass the LocID and dbpath if a special case. A new
        Enceladus will be initialised with the parameters in the database."""

        db = sqlite3.connect(dbpath)
        cursor = db.cursor()
        # select the row for your methanogen
        cursor.execute("SELECT * FROM Venus WHERE LocID = ?", (LocID,))
        params = cursor.fetchone()
        CompID = params[4]
        db.close()

        E = cls(LocID, H2ppm=25, setupcomp=False)

        c, pH = VenusDrop.db_to_comp(CompID, E.env)

        # print(c)
        E.composition = c
        E.initial_conditions(25, setupcomp=False)

        return E

    @staticmethod
    def db_to_comp(CompID, envir, dbpath=std_dbpath):
        db = sqlite3.connect(dbpath)
        cursor = db.cursor()
        cursor.execute("SELECT * FROM composition WHERE CompID = ?", (CompID,))

        params = cursor.fetchone()
        Comp = ast.literal_eval(params[1])
        pH = params[-1]

        db.close()

        # add in CHNOPS, in elemental form?

        H2 = reaction.reagent('H2(aq)', self.env, phase='aq' ,conc=Comp['H2(aq)'], activity=Comp['H2(aq)'])
        SO4 = reaction.reagent('SO4--(aq)', self.env, phase='aq', charge=-2, conc=Comp['SO4--(aq)'], activity=Comp['SO4--(aq)'])
        H = reaction.reagent('H+(aq)', self.env, phase='aq', charge=1, conc=Comp['H+(aq)'], activity=Comp['H+(aq)'])
        HS = reaction.reagent('H2(aq)', self.env, phase='aq' , charge=-1, conc=mComp['H2(aq)'], activity=Comp['H2(aq)'])

        H2O = reaction.reagent('H2O(l)', self.env, phase='l', conc=55.5,
          phase_ss=True, activity=1)


        composition = {H2.name:H2, SO4.name:SO4,
          H.name:H, H2O.name:H2O, HS.name:HS}

        return composition, pH

    def to_db(self, updateComp=True, dbpath=std_dbpath):
        """Send an Enceladus to be stored in the database"""
        pH = -math.log10(self.composition['H+'].conc)
        if updateComp:
            self.Comp_to_db(CompID, pH)

        db = sqlite3.connect(dbpath)
        cursor = db.cursor()

        # fill in the methanogen database
        try:
            cursor.execute(" INSERT INTO Venus(LocID, name, pressure, temperature, CompID, volume) VALUES(?,?,?,?,?,?)", ('Tester', self.name, self.env.P, self.env.T, self.CompID, self.env.V))

            #if that goes through, replace tester with a generated OrgID.
            cursor.execute('SELECT rowid FROM Venus WHERE LocID = ?', ('Tester',))
            entryno = cursor.fetchone()[0]
            self.LocID = self.name+str(entryno)+date.today().strftime("_%d%m%y") # ecosystem style
            print(self.LocID)
            cursor.execute('UPDATE Venus SET LocID = ? WHERE LocID = ?', (self.LocID, 'Tester'))

            db.commit()

        except sqlite3.IntegrityError as e:
            print(str(e))
            if 'column LocID' in str(e):#.startswith('column OrgID'):
                print("\n LocID already in use. Maybe Tester hasn't been deleted somewhere \n")
                raise e
            elif 'columns' or 'UNIQUE' in str(e):#.startswith('columns'):
                print('\n A Venus with these parameters already exists in the database \n')

                cursor.execute("SELECT LocID FROM Venus WHERE name = ? AND pressure = ? AND temperature = ? AND CompID = ? AND volume = ?", (self.name, self.env.P, self.env.T, self.CompID, self.env.V))

                ID = cursor.fetchone()
                print("\n Setting this Venus' ID to " + ID[0] + "\n")
                self.LocID = ID[0] #set it here
            else:
                raise e
        finally:
            #ensure we safely close the database even if there is an error.
            db.close()
            print(self.LocID)

        #now its in, send the reactor to the locale database
        self.update_locale_db(dbpath)


    def workoutID(self, dbpath=std_dbpath):
        #first we need the compID for the present compostion
        pH = -math.log10(self.composition['H+'].conc)
        self.Comp_to_db(pH)

        db = sqlite3.connect(dbpath)
        cursor = db.cursor()

        # search in the database for Encealdus' unique identifiers

        try:
            cursor.execute('SELECT LocID FROM Venus WHERE name = ? AND pressure = ? AND temperature = ? AND CompID = ? AND volume = ? AND pH = ? AND depth = ?', (self.name, self.env.P, self.env.T, self.CompID, self.env.V, -math.log10(self.composition['H+'].conc), self.depth))

            op = cursor.fetchone()
            if op is not None:
                print("\n Venus found! It's LocID is: "+op[0]+" \n")
                print(op[0])
                self.LocID = str(op[0])
            else:
                raise Exception('No Venus found')
        except Exception as e:
            print(e.args)
            print('adding to database')
            self.to_db(updateComp=False, dbpath=dbpath)
        finally:
            db.close()


    def initial_conditions(self, H2ppm, HSact, setupcomp=True):

        if setupcomp:
            # From Waite:2017, estimate the molality of CO2 from the pH
            mol_H2 = VenusDrop.getgasconc('H2', self.env.P*H2ppm/(1e6), self.env.T, P_bar=self.env.P)

            mol_SO4 = 0.5 # assume fully dissoiated at pH 0 (from H2SO4)
            mol_H = 1.0 # pH zero
            mol_HS = 0.0 #? don't know, assume biogenic only.

            SO4 = reaction.reagent('SO4--', self.env, phase='aq', charge=-2, conc=mol_SO4, activity=mol_SO4)
            H2 = reaction.reagent('H2(aq)', self.env, phase='aq' ,conc=mol_H2, activity=mol_H2)

            H = reaction.reagent('H+', self.env, phase='aq', charge=1, conc=mol_H, activity=mol_H)
            HS = reaction.reagent('HS-', self.env, phase='aq' , charge=-1, conc=HSact, activity=HSact)

            H2O = reaction.reagent('H2O(l)', self.env, phase='l', conc=55.5,
              phase_ss=True, activity=1)


            self.composition = {H2.name:H2, SO4.name:SO4,
              H.name:H, H2O.name:H2O, HS.name:HS}

            # add in CHNOPS, in elemental form


        # el = reaction.reagent('e-', self.env, charge=-1)
        # # overall
        r = {self.composition['H2(aq)']:4, self.composition['SO4--']:1, self.composition['H+']:1}
        p = {self.composition['HS-']:1, self.composition['H2O(l)']:4}
        thermaloa = reaction.reaction(r,p,self.env)


        # # redox
        # fr = {self.composition['CO2(aq)']:1, self.composition['H+']:8, el:8}
        # fp = {self.composition['CH4(g)']:1, self.composition['H2O(l)']:2}
        # fwd = redox_half(fr, fp, self.env, 8, -0.244)
        #
        # rr = {self.composition['H+']:8, el:8}
        # rp = {self.composition['H2(aq)']:4}
        #
        # # E ref: Karp, Gerald, Cell and Molecular Biology, 5th Ed., Wiley, 2008
        # rvs = redox_half(rr, rp, self.env, 8, -0.421)
        #
        # # redox cell reaction is thermal overall
        # rx = redox(fwd, rvs, thermaloa, 8, self.env)

        self.reactionlist = {'SulfateReduction':thermaloa}

    @staticmethod
    def solubilityEQparams(A1, A2, A3, B1, B2, B3, T, S):
        lnc = A1+(A2*100/T)+(A3*math.log(T/100)) + (S*(B1+(B2*T/100)+(B3*((T/100)**2))))
        return math.exp(lnc)

    @staticmethod
    def solubilityEQ(gas_name, T, S):
        soldict = gases[gas_name]['SolubilityCons']
        return VenusDrop.solubilityEQparams(soldict[0], soldict[1], soldict[2], soldict[3], soldict[4], soldict[5], T, S)

    @staticmethod
    def getgasconc(gas_name, P_gas, T, P_bar=101325, S=0):
        """P_gas ---partial pressure of gas in Pa"""

        Beta = VenusDrop.solubilityEQ(gas_name, T, S) #will end up in mL / mL atm ? or should
        if gases[gas_name]['SolubilityForm']=='K':
            #with bar on the end is risky
            Beta = Beta * R*298.15*P_bar/101325


        #calculate the volume of dissolved gas found per liter of water:
        C_gas = (P_gas/101325)*Beta
        #calculate no of moles per L of gas from ideal gas law:
        M_gas = (P_gas/101325)/(R*T)
        conc_gas = C_gas*M_gas
        return(conc_gas)
