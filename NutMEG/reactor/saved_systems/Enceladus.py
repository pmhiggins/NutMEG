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

    """
    The chemically active evirionment of the Enceladean Ocean.
    It has a few unique attributes, and at the moment does not take
    kwargs for reactors. This may be added in the future when the
    model is a bit more flexible.

    This is a saved instance of reactor, which contains the
    methanogenesis and environmental data for the moon. In the future
    this will be fleshed out with other constituents and maybe
    even other reactions/flows.

    Attributes
    ----------
    ocean_pH : float
        pH of the ocean at 273.15 K
    mixingratios : dict, optional
        molar mixing ratio of the plume (wrt H2O) to use as ufloats. Default is
        the mixing ratios listed in Waite et al. 2017.
    CO2origin : str
        By which methods to calculate the activity of CO2. Options include
        'tigerstripe', 'pH', 'HTHeating', 'HTHeating20', 'HTHeatingSalts'.
        HTHeatingSalts is recommended for the best results.
    tigerstripeT : ufloat, optional
        Temperature of the tiger stripes for estimating CO2 activity. Only
        required if you plan to use this technique
        (e.g. when CO2origin = 'tigerstripe')
    nominals : boolean, optional
        Whether to return values as their nominal case (True) or as ufloats
        (False). Using ufloats can cause problems when working out reaction
        quotients.
    Pconc : ufloat, optional
        Dissolved concentration of phosphorus in the ocean.
    workoutID : boolean, optional
        Whether to make an entry into the Enceladus table of the NutMEG
        database.

    """

    volume = 7.54e15 #m^3, from radius = 250km, depth = 10km
    env = environment(T=300., P=70e5, V=volume)
    depth=8.5


    def __init__(self, name, pH=8.5, depth=8.5, T=273.15,
      mixingratios={},
      CO2origin='pH',
      tigerstripeT=uf(197,20),
      Pconc=uf(1e-6,0), workoutID=False,
      nominals=False,
      saltlevel='nom',
      **kwargs):
        self.depth=depth
        self.env.T=T
        self.env.P=(10.+(7.*self.depth))*100000 #80 at bottom, 10 at top, linear in between.
        self.ocean_pH = pH # the pH of the 'wider ocean' e.g. at 273 K.
        self.pH = pH # the pH at this section of the ocean, not nec. same as above.
        self.mixingratios=Waite2017ratios
        self.mixingratios.update(mixingratios)
        self.tigerstripeT = tigerstripeT
        self.nominals=nominals

        mol_CO2 = 0.
        aH2O=1.
        oceanvals=False

        if CO2origin=='tigerstripe':
            mol_CO2 = self.get_tigerstripe_CO2()
            self.initial_conditions(self.pH, mol_CO2, Pconc, H2Oact=aH2O, oceanvals=oceanvals)
        elif CO2origin=='pH':
            mol_CO2 = self.get_CO2_from_pH()
            self.initial_conditions(self.pH, mol_CO2, Pconc, H2Oact=aH2O, oceanvals=oceanvals)
        elif CO2origin=='HTHeating':
            mol_CO2, aH2O = Enceladus.get_CO2_from_HTHeating(T=self.env.T, pH=self.ocean_pH, nominals=self.nominals)
            oceanvals=True
            self.initial_conditions(self.pH, mol_CO2, Pconc, H2Oact=aH2O, oceanvals=oceanvals)
        elif CO2origin=='HTHeating20':
            mol_CO2, aH2O = Enceladus.get_CO2_from_HTHeating(T=self.env.T, pH=self.ocean_pH, nominals=self.nominals, CO2unc=0.2)
            oceanvals=True
            self.initial_conditions(self.pH, mol_CO2, Pconc, H2Oact=aH2O, oceanvals=oceanvals)
        elif CO2origin=='HTHeatingSalts':
            mol_CO2lst, aH2Olst = Enceladus.get_CO2_from_HTHeating(T=self.env.T, pH=self.ocean_pH, nominals=self.nominals, CO2unc=0., salts=True)
            # mol_CO2, aH2O = mol_CO2lst[0], aH2Olst[0]
            oceanvals=True
            if saltlevel == 'nom':
                self.initial_conditions(self.pH, mol_CO2lst[0], Pconc, H2Oact=aH2Olst[0], oceanvals=oceanvals)
            elif saltlevel == 'high':
                self.initial_conditions(self.pH, mol_CO2lst[1], Pconc, H2Oact=aH2Olst[1], oceanvals=oceanvals)
            elif saltlevel == 'low':
                self.initial_conditions(self.pH, mol_CO2lst[2], Pconc, H2Oact=aH2Olst[2], oceanvals=oceanvals)

        reactor.__init__(self, name, env=self.env,
          reactionlist=self.reactionlist,
          composition=self.composition,
          workoutID=workoutID, **kwargs)

    def get_tigerstripe_CO2(self, logform=False):
        """
        Return the CO2 activity as estimated by the tiger stripe temperature.
        (Glein and Waite 2020)
        """
        if logform:
            return umath.log10(self.mixingratios['CO2']/100)+9.429-(2574/self.tigerstripeT)
        else:
            return umath.pow(10, (umath.log10(self.mixingratios['CO2']/100)+9.429-(2574/self.tigerstripeT)))

    def get_CO2_from_pH(self):
        """
        Return the CO2 activity as estimated from the pH in (Waite et al. 2017).
        Note only really valid at 273 K.
        """
        absCO2 = 10**(-0.1213*(self.pH*self.pH) + (0.9832*self.pH) -3.1741)
        unc = 0.2*absCO2
        if self.nominals:
            return uf(absCO2, unc).n
        else:
            return uf(absCO2, unc)

    @staticmethod
    def interpo_CO2_H2O(Temp=273.15, oceanpH=8.0, salt='nominalCO2'):
        """
        Return the activity of CO2 and H2O at a given temperature and
        wider ocean pH (273 K) by interpolating the values from the
        carbonate speciation model.
        """
        fn_preamble = os.path.dirname(__file__)+'/../../data/Enceladus/'

        aH2O = np.load(fn_preamble+salt+'/aH2O.npy')
        aCO2 = np.load(fn_preamble+salt+'/aCO2.npy')
        # pHHT = np.load('HTHeatingdata/nominalCO2/pHHT.npy')

        pHfloats = np.linspace(7.,12., num=11)
        Tfloats = np.linspace(273.15, 473.15, num=21)

        fCO2 = interpolate.interp2d(Tfloats,pHfloats,aCO2,kind='cubic')

        fH2O = interpolate.interp2d(Tfloats,pHfloats,aH2O,kind='cubic')

        return fCO2(Temp, oceanpH)[0], fH2O(Temp, oceanpH)[0]

    @staticmethod
    def get_CO2_from_HTHeating(T, pH, nominals=False, salts=False, CO2unc=0.):
        """
        Return the CO2 activity and H2O activity as ufloats from the carbonate
        speciation model, or as a list with upper and lower bounds
        (when salts==True)

        Attributes
        ----------
        T : float
            Temperature of your position in the ocean
        pH : float
            pH that the ocean would be at 273 K.
        salts : boolean, optional
            if True, use the carbonate speciation model and return
            nominal, upper and lower bounds for both.
        """
        aCO2, aH2O = Enceladus.interpo_CO2_H2O(T, pH)
        if aCO2<0:
            print('WARNING: negative CO2 at pH ', pH, 'Temperature', T, 'setting to 1e10 M')
            if nominals:
                return 1e-10, aH2O
            elif salts:
                aCO2_hs, aH2O_hs = Enceladus.interpo_CO2_H2O(T, pH, salt='highsalt')
                aCO2_ls, aH2O_ls = Enceladus.interpo_CO2_H2O(T, pH, salt='lowsalt')
                return [1e-10, 1e-10, 1e-10], [aH2O, aH2O_hs, aH2O_ls]
            else:
                return uf(1e-10, 2e-11), uf(aH2O,0)
        if nominals:
            return aCO2, aH2O
        elif salts:
            aCO2_hs, aH2O_hs = Enceladus.interpo_CO2_H2O(T, pH, salt='highsalt')
            aCO2_ls, aH2O_ls = Enceladus.interpo_CO2_H2O(T, pH, salt='lowsalt')
            return [aCO2, aCO2_hs, aCO2_ls], [aH2O, aH2O_hs, aH2O_ls]
        else:
            return uf(aCO2, CO2unc*aCO2), uf(aH2O,0)



    @staticmethod
    def getHsuT(pH):
        """
        Get the temperature at the rock-water interface predicted by
        Hsu et al 2015 if pH were constant.
        """
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
        """ Use simple exponential model to estimate T with depth."""
        return (T_max-T_min)*math.exp(-1.*height)+T_min

    def calc_mol_CH4(self, mol_CO2_oc):
        """
        Return the activity of CH4 from the activity of CO2 and mixing ratios.
        """
        return (self.mixingratios['CH4']/self.mixingratios['CO2'])*mol_CO2_oc

    def calc_mol_H2(self, mol_CO2_oc):
        """
        Return the activity of H2 from the activity of CO2 and mixing ratios.
        """
        return (self.mixingratios['H2']/self.mixingratios['CO2'])*mol_CO2_oc


    def initial_conditions(self, pH, mol_CO2, Pconc, H2Oact=1.0, oceanvals=False, mol_CO2_oc=None):
        """
        Set up the initial conditions of the ocean in the configuration
        defined in the initialisation.

        Attributes
        ----------
        pH : float
            wider ocean (273 K) pH
        mol_CO2 : float
            Calculated activity of CO2 at this ocean location.
        Pconc : ufloat
            Concentration of phosphorus
        H2Oact : float, optional
            water activity, default 1.0.
        oceanvals : boolean, optional
            whether computing the wider ocean CO2 activity is needed.
        mol_CO2_oc : float, optional
            wider ocean CO2 activity. Default is None (then calculated here)
        """
        if mol_CO2_oc == None:
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
        CH4aq = reaction.reagent('Methane(aq)', self.env, phase='g', conc=mol_CH4,
          activity=mol_CH4)
        H2O = reaction.reagent('H2O(l)', self.env, phase='l', conc=uf(55.5, 0), activity=uf(H2Oact,0))
        el = reaction.reagent('e-', self.env, charge=-1)
        H = reaction.reagent('H+', self.env, charge=1, conc=mol_H,
          phase='aq', activity=mol_H)

        self.composition = {CO2.name:CO2, H2aq.name:H2aq,
          CH4aq.name:CH4aq, H2O.name:H2O, H.name:H}

        # overall
        r = {self.composition['CO2(aq)']:1, self.composition['H2(aq)']:4}
        p = {self.composition['Methane(aq)']:1, self.composition['H2O(l)']:2}
        thermaloa = reaction.reaction(r,p,self.env)

        # redox
        fr = {self.composition['CO2(aq)']:1, self.composition['H+']:8, el:8}
        fp = {self.composition['Methane(aq)']:1, self.composition['H2O(l)']:2}
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
        # P and S we don't actually know.
        self.composition['P(aq)'] = reaction.reagent('P(aq)', self.env, phase='aq', conc=Pconc,
          activity=Pconc, thermo=False)
        self.composition['H2S(aq)'] = reaction.reagent('H2S(aq)', self.env, phase='aq', conc=mol_H2S,
          activity=mol_H2S, thermo=False)
