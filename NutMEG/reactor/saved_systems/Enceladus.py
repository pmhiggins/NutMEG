import sys, os
sys.path.append(os.path.dirname(__file__))

import math
from NutMEG.environment import *
from NutMEG.reactor import reactor
import NutMEG.reaction as reaction
from itertools import chain

from uncertainties import ufloat as uf
import uncertainties.umath as umath
import numpy as np
from scipy import interpolate



Waite2017ratios = {'CO2': uf(0.55, 0.25), 'CH4': uf(0.2, 0.1), 'NH3': uf(0.85, 0.45), 'H2': uf(0.9,0.5), 'H2S':uf(0.0021,0.001)}

class Enceladus(reactor):

    """The chemically active evirionment of the Enceladean Ocean.

    This is a saved instance of reative_system, which contains the
    methanogenesis and envrironmental data for the moon. In the full
    model this will be fleshed out with other constituents and maybe
    even other reactions/flows.
    """

    volume = 7.54e15 #m^3, from radius = 250km, depth = 10km
    env = environment(T=300., P=70e5, V=volume)
    depth=8.5


    def __init__(self, name, pH=8.5, depth=8.5, T=273.15,
      mixingratios=Waite2017ratios,
      CO2origin='pH',
      tigerstripeT=uf(197,20),
      Pconc=uf(1e-6,0), workoutID=False,
      nominals=False):
        self.depth=depth
        self.env.T=T
        self.env.P=(10.+(7.*self.depth))*100000 #80 at bottom, 10 at top, linear in between.
        self.pH = pH
        self.mixingratios=mixingratios
        self.tigerstripeT = tigerstripeT
        self.nominals=nominals

        mol_CO2 = 0.
        aH2O=1.
        oceanvals=False

        if CO2origin=='tigerstripe':
            mol_CO2 = self.get_tigerstripe_CO2()
        elif CO2origin=='pH':
            mol_CO2 = self.get_CO2_from_pH()
        elif CO2origin=='HTHeating':
            mol_CO2, aH2O = Enceladus.get_CO2_from_HTHeating(T=self.env.T, pH = self.pH, nominals=self.nominals)
            oceanvals=True
        elif CO2origin=='HTHeating20':
            mol_CO2, aH2O = Enceladus.get_CO2_from_HTHeating(T=self.env.T, pH = self.pH, nominals=self.nominals, CO2unc=0.2)
            oceanvals=True

        self.initial_conditions(self.pH, mol_CO2, Pconc, H2Oact=aH2O, oceanvals=oceanvals)

        reactor.__init__(self, name, env=self.env,
          reactionlist=self.reactionlist,
          composition=self.composition,
          workoutID=workoutID)

    def get_tigerstripe_CO2(self, logform=False):
        if logform:
            return umath.log10(self.mixingratios['CO2']/100)+9.429-(2574/self.tigerstripeT)
        else:
            return umath.pow(10, (umath.log10(self.mixingratios['CO2']/100)+9.429-(2574/TigerstripeT)))

    def get_CO2_from_pH(self):
        absCO2 = 10**(-0.1213*(self.pH*self.pH) + (0.9832*self.pH) -3.1741)
        unc = 0.2*absCO2
        if self.nominals:
            return uf(absCO2, unc).n
        else:
            return uf(absCO2, unc)

    @staticmethod
    def interpo_CO2_H2O(Temp=273.15, oceanpH=8.0):

        aH2O = np.load(os.path.dirname(__file__)+'/../../data/Enceladus/aH2O.npy')
        aCO2 = np.load(os.path.dirname(__file__)+'/../../data/Enceladus/aCO2.npy')
        # pHHT = np.load('HTHeatingdata/nominalCO2/pHHT.npy')

        pHfloats = np.linspace(7.,12., num=11)
        Tfloats = np.linspace(273.15, 473.15, num=21)

        fCO2 = interpolate.interp2d(Tfloats,pHfloats,aCO2,kind='cubic')

        fH2O = interpolate.interp2d(Tfloats,pHfloats,aH2O,kind='cubic')

        return fCO2(Temp, oceanpH)[0], fH2O(Temp, oceanpH)[0]

    @staticmethod
    def get_CO2_from_HTHeating(T, pH, nominals=False, CO2unc=0.):
        aCO2, aH2O = Enceladus.interpo_CO2_H2O(T, pH)
        if nominals:
            return aCO2, aH2O
        else:
            return uf(aCO2, CO2unc*aCO2), uf(aH2O,0)



    @staticmethod
    def getHsuT(pH):
        """Get the temperature at the rock-water interface predicted by
        Hsu et al 2015 if pH were constant"""
        if pH==8:
            return 230.0+273.15
        if pH>8 and pH<9:
            return (230-((pH-8)*(230-147)))+273.15
        if pH==9:
            return 147+273.15
        if pH>9 and pH<10:
            return (147-((pH-9)*(147-97)))+273.15
        if pH==10:
            return 97+273.15
        if pH>10 and pH<11:
            return (97-((pH-10)*(97-90)))+273.15
        if pH>=11:
            return 90+273.15

    def get_exp_T(self, T_max, T_min, height):
        return (T_max-T_min)*math.exp(-1.*height)+T_min



    def initial_conditions(self, pH, mol_CO2, Pconc, H2Oact=1.0, oceanvals=False):

        mol_CO2_oc = mol_CO2
        if oceanvals:
            #get wider ocean CO_2
            mol_CO2_oc = Enceladus.get_CO2_from_HTHeating(T=273.15, pH=pH)[0]

        mol_CH4 = (self.mixingratios['CH4']/self.mixingratios['CO2'])*mol_CO2_oc#2.75*mol_CO2 #0.34*mol_CO2
        mol_H2 = (self.mixingratios['H2']/self.mixingratios['CO2'])*mol_CO2_oc#11*mol_CO2 #103*mol_CO2
        mol_NH3 = (self.mixingratios['NH3']/self.mixingratios['CO2'])*mol_CO2_oc
        mol_H2S = (self.mixingratios['H2S']/self.mixingratios['CO2'])*mol_CO2_oc
        # ^ from the ratios in the plumes, hence are upper limits
        mol_H=uf(10**(-pH),0)
        if self.nominals:
            mol_CH4 = mol_CH4.n
            mol_H2 = mol_H2.n
            mol_NH3 = mol_NH3.n
            mol_H2S = mol_H2S.n
            mol_H=uf(10**(-pH),0).n

        # reagents
        CO2 = reaction.reagent('CO2(aq)', self.env, phase='g', conc=mol_CO2,
          activity=mol_CO2)
        H2aq = reaction.reagent('H2(aq)', self.env, phase='aq', conc=mol_H2,
          activity=mol_H2)
        CH4aq = reaction.reagent('CH4(g)', self.env, phase='g', conc=mol_CH4,
          activity=mol_CH4)
        H2O = reaction.reagent('H2O(l)', self.env, phase='l', conc=uf(55.5, 0), activity=uf(1,0))
        el = reaction.reagent('e-', self.env, charge=-1)
        H = reaction.reagent('H+', self.env, charge=1, conc=mol_H,
          phase='aq', activity=mol_H)

        self.composition = {CO2.name:CO2, H2aq.name:H2aq,
          CH4aq.name:CH4aq, H2O.name:H2O, H.name:H}

        # overall
        r = {self.composition['CO2(aq)']:1, self.composition['H2(aq)']:4}
        p = {self.composition['CH4(g)']:1, self.composition['H2O(l)']:2}
        thermaloa = reaction.reaction(r,p,self.env)

        # redox
        fr = {self.composition['CO2(aq)']:1, self.composition['H+']:8, el:8}
        fp = {self.composition['CH4(g)']:1, self.composition['H2O(l)']:2}
        fwd = reaction.redox_half(fr, fp, self.env, 8, -0.244)

        rr = {self.composition['H+']:8, el:8}
        rp = {self.composition['H2(aq)']:4}

        # E ref: Karp, Gerald, Cell and Molecular Biology, 5th Ed., Wiley, 2008
        rvs = reaction.redox_half(rr, rp, self.env, 8, -0.421)

        # redox cell reaction is thermal overall
        rx = reaction.redox(fwd, rvs, thermaloa, 8, self.env)

        self.reactionlist = {thermaloa.equation:{type(thermaloa):thermaloa,
          type(rx):rx}}

        # add in CHNOPS, in elemental form for now
        # we already have C in the form of CO2 and CH4
        self.composition['NH3(aq)'] = reaction.reagent('NH3(aq)', self.env, phase='aq', conc=mol_NH3,
          activity=mol_NH3) # from Glein, Baross, Waite 2015
        # plenty of O in water lol
        # H is also in water, and the dissolved H2
        # P and S we dont actually know, so let's arbitrarily say 1 micromole
        self.composition['P(aq)'] = reaction.reagent('P(aq)', self.env, phase='aq', conc=Pconc,
          activity=Pconc, thermo=False)
        self.composition['H2S(aq)'] = reaction.reagent('H2S(aq)', self.env, phase='aq', conc=mol_H2S,
          activity=mol_H2S, thermo=False)






""" These are probably redundant but I'll hold onto them for now.

    def perform_methanogenesis(self, molesCO2, number=1):
        for c in self.reactionlist['Methanogenesis(thermo)'].reactants:
            # Find total number of moles in system, then remove the amount
            # that has been reacted away or formed.
            c.molal = (((c.activity*1000.0*self.volume)
              - (c.molar_ratio*molesCO2*number))/(1000.0*self.volume))
            c.activity = c.molal
        for c in self.reactionlist['Methanogenesis(thermo)'].products:
            # Find total number of moles in system, then remove the amount
            # that has been reacted away or formed.
            c.molal = (((c.activity*1000.0*self.volume)
              + (c.molar_ratio*molesCO2*number))/(1000.0*self.volume))
            c.activity = c.molal

    def perform_methanogenesis_redox(self, molesCO2, number=1):
        for c in chain(self.reactionlist['Methanogenesis(redox)'].forward.reactants,
          self.reactionlist['Methanogenesis(redox)'].reverse.products):
            # Find total number of moles in system, then remove the amount
            # that has been reacted away or formed.
            c.molal = (((c.activity*1000.0*self.volume)
              - (c.molar_ratio*molesCO2*number))/(1000.0*self.volume))
            c.activity = c.molal
        for c in chain(self.reactionlist['Methanogenesis(redox)'].reverse.reactants,
          self.reactionlist['Methanogenesis(redox)'].forward.products):
            # Find total number of moles in system, then remove the amount
            # that has been reacted away or formed.
            c.molal = (((c.activity*1000.0*self.volume)
              + (c.molar_ratio*molesCO2*number))/(1000.0*self.volume))
            c.activity = c.molal
"""
