import sys, os
sys.path.append('..')
import math

# NutMEG imports
import NutMEG
from NutMEG.reactor import reactor
import NutMEG.reaction as reaction


R=0.08206 # L atm / mol K, universal gas constant.

# some solubility data for various gases.
gases = {
  'H2(aq)': {
    'SolubilityCons':[-39.9611, 53.9381, 16.3135, -0.036249,
      0.017565,-0.0023010], 'SolubilityForm':'Beta'},
  'CO2(aq)' : {
    'SolubilityCons':[-58.0931, 90.5069, 22.2940, 0.027766,
      -0.025888,0.0050578], 'SolubilityForm':'K'},
  'O2(aq)' : {
    'SolubilityCons':[-58.3877, 85.8079, 23.8439, -0.034892,
      0.015568,-0.0019387], 'SolubilityForm':'Beta'},
  'CH4(g)' : {
    'SolubilityCons':[-68.8862, 101.4956, 28.7314, -0.076146,
      0.043970,-0.006872], 'SolubilityForm':'Beta'}}


class VenusDrop(reactor):

    """The chemically active evironment of a simple Venusian cloud droplet.

    This is an instance of reactor, with added methods to estimate the
    composition from cloud ppm values and setup specific possible venusian
    metabolisms.
    """

    def __init__(self, H2ppm=25, HSact=1e-10, T=298, P=100000, V=1.13e-16,
      setupcomp=True, workoutID=True):

        self.env = NutMEG.environment(T=T, P=P, V=V)
        self.name='VenusDrop'

        if setupcomp:
            self.composition={}
            self.reactionlist={}
            self.initial_conditions(H2ppm, HSact=HSact,
              setupcomp=setupcomp)
            reactor.__init__(self, self.name, env=self.env,
              reactionlist=self.reactionlist,
              composition=self.composition,
              pH=self.pH, workoutID=workoutID)
        else:
            # the conditions will be set up in the test code,
            # so initiailise a default reactor
            reactor.__init__(self, self.name, env=self.env,
              workoutID=workoutID)



    def initial_conditions(self, H2ppm, HSact=1e-10, setupcomp=True, addreaction=True):
        """ Set up a basic configuration, with pH 0 sulfuric acid
        droplet and dissolved H2.

        H2ppm : float
            H2 parts per million in the atmosphere around the droplet.
        HSact : float
            HS activity to use. Default is 1e-10
        setupcomp : bool
            Whether to set up the standard composition. Default True
        addreaction : bool
            Whether to include sulfate reduction in the reactionlist
        """

        if setupcomp:

            mol_H2 = VenusDrop.getgasconc('H2(aq)', self.env.P*H2ppm/(1e6),
              self.env.T, P_bar=self.env.P)
            mol_SO4 = 0.5 # assume fully dissoiated at pH 0 (from H2SO4)
            mol_H = 1.0 # pH zero
            mol_HS = HSact #? don't know, let's assume biogenic only.


            SO4 = reaction.reagent('SO4--', self.env, phase='aq', charge=-2,
              conc=mol_SO4, activity=mol_SO4)
            H2 = reaction.reagent('H2(aq)', self.env, phase='aq',
              conc=mol_H2, activity=mol_H2)

            H = reaction.reagent('H+', self.env, phase='aq', charge=1,
              conc=mol_H, activity=mol_H)
            HS = reaction.reagent('HS-', self.env, phase='aq', charge=-1,
              conc=HSact, activity=HSact)

            H2O = reaction.reagent('H2O(l)', self.env, phase='l', conc=55.5,
              phase_ss=True, activity=1)

            # add these reagents to the composition of the VenusDrop
            self.composition.update({H2.name:H2, SO4.name:SO4,
              H.name:H, H2O.name:H2O, HS.name:HS})


        self.pH = 10**(self.composition['H+'].activity)

        if addreaction:
            # put together an overall reaction fro sulfate reduction
            r = {self.composition['H2(aq)']:4, self.composition['SO4--']:1,
              self.composition['H+']:1}
            p = {self.composition['HS-']:1, self.composition['H2O(l)']:4}
            thermaloa = reaction.reaction(r,p,self.env)

            # add this reaction to the VenusDrop reactor.
            self.add_reaction(thermaloa, overwrite=True)

    def update_conditions(self, H2ppm, HSact=1e-10):
        """Calculate dissolved concentratoins in the current reactor
        environment for a standard VenusDrop

        H2ppm : float
            Atmospheric ppm for H2
        """

        molH2 = VenusDrop.getgasconc('H2(aq)', self.env.P*H2ppm/(1e6),
          self.env.T, P_bar=self.env.P)

        self.composition['H2(aq)'].activity = molH2
        self.composition['H2(aq)'].conc = molH2
        self.composition['HS-'].activity = HSact
        self.composition['HS-'].conc = HSact


    def update_reagent(self, name, ppm):
        """ Update the concentration of a specific reagent based
        on its atmospheric ppm

        name : str
            Name of the reagent to update
        ppm : float
            atmospheric concentration to use
        """

        mol_r = VenusDrop.getgasconc(name, self.env.P*ppm/(1e6),
            self.env.T, P_bar=self.env.P)

        r = reaction.reagent(name, self.env, phase='aq', charge=0,
              conc=mol_r, activity=mol_r)

        self.composition.update({r.name:r})



    @staticmethod
    def solubilityEQparams(A1, A2, A3, B1, B2, B3, T, S):
        lnc = (A1+(A2*100/T)+(A3*math.log(T/100)) +
          (S*(B1+(B2*T/100)+(B3*((T/100)**2)))))
        return math.exp(lnc)

    @staticmethod
    def solubilityEQ(gas_name, T, S):
        soldict = gases[gas_name]['SolubilityCons']
        return VenusDrop.solubilityEQparams(soldict[0], soldict[1],
          soldict[2], soldict[3], soldict[4], soldict[5], T, S)

    @staticmethod
    def getgasconc(gas_name, P_gas, T, P_bar=101325, S=0):
        """From a partial pressure, get the dissolved gas
        concentration.

        gas_name : str
            Name of reagent
        P_gas : float
            partial pressure of gas in Pa
        T : float
            local temperature in K
        P_bar : float
            Local barometric pressure in Pa
        S : float
            Salinity of the droplet in parts per thousand

        """

        Beta = VenusDrop.solubilityEQ(gas_name, T, S)
        # units mL / mL atm

        if gases[gas_name]['SolubilityForm']=='K':
            Beta = Beta * R*298.15*P_bar/101325

        #calculate the volume of dissolved gas found per liter of water:
        C_gas = (P_gas/101325)*Beta
        #calculate no of moles per L of gas from ideal gas law:
        M_gas = (P_gas/101325)/(R*T)
        conc_gas = C_gas*M_gas
        return conc_gas
