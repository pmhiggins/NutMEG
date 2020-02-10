"""

This is the maintainer class, which computes and calculates the maintenance
requirements in the form of powers for a given organism.

IMPORTANT: All values are PER CELL, so for a Horde, when you use these
attributes you will likely need to scale them up.

"""
from .adaptations.pHadaptations import pHadaptations
from .adaptations.Tadaptations import Tadaptations
import ast

import NutMEG.util.NutMEGparams as nmp
from NutMEG.util.loggersetup import loggersetup as logset
logger = logset.get_logger(__name__, filelevel=nmp.filelevel, printlevel=nmp.printlevel)

class maintainer:

    P_loss = 0. # Instantaneous fractional power loss due to all maintenance.
    net_dict = {} # loss mechanisms, each in W
    frac_dict = {}  # loss mechanisms as fractions of supply

    P_store = 0. # Fractional power stored away for future use.
    #P_growth = 0. # Total power left over --- put into growth.


    def __init__(self, host,
      net_dict={},
      Basal=0.0,
      supply=1.0,
      Tdef='None',
      pHdef='None'):

        self.host = host # this is a reference to the host organism for
          # e.g. if the environment changes.
        self.net_dict.update(net_dict)
        self.net_dict['Basal']=Basal
        self.Tdef = Tdef
        self.get_P_T()
        self.pHdef = pHdef
        self.get_P_pH()
        self.update_frac_dict(supply)


    def add_mechanism(self, name, val):
        """Add a mechanism to the dictionary of coping mechanisms."""
        self.net_dict[name] = val


    def update_frac_dict(self, P_supply):#, hosthorde=True):
        """Update the power loss dictionary by expressing it as a fraction
        of power supply.
        """
        if P_supply != 0:
            for key, value in self.net_dict.items():
                self.frac_dict[key] = (value/P_supply)
        else:
            for key, value in self.net_dict.items():
                self.frac_dict[key] = (0)
            # hacky, stops ValueErrors later in the step
            self.host.E_growth=0


    def update_P_loss(self, P_supply):
        """Use the power loss dictionary to get the current total
        power lost to maintenance.
        """

        self.update_frac_dict(P_supply)
        loss = 0.

        for val in self.frac_dict.values():
            loss += (val)#*0.001*random.randrange(500,1500))
        if self.host.locale.env.T > (400.): #proteins definitely break down
            self.P_loss = 1.0
        elif loss<=1.0:
            self.P_loss = loss
        else:
            self.P_loss = 1.0



    def get_P_T(self):
        """ Calculate the power cost related to temperature and update
        the net_dict. Make sure you have set Tdef to the defence that
        you want to compute!
        """
        T_ad_calc = Tadaptations(self.host)
        if self.Tdef=='Lever10pc':
            # Use the Lever calculation with replacement at 10% [QE]
            self.net_dict['T'] = T_ad_calc.getLeverME()
        elif self.Tdef=='Lever2pc':
            # Use the Lever calculation with replacement at 2% [QE]
            self.net_dict['T'] =  T_ad_calc.getLeverME(cutoff_pc=2)
        elif self.Tdef=='Tijhuis':
            # use the Tijhuis calculation [QE]
            self.net_dict['T'] = T_ad_calc.getTijuisME()
        elif self.Tdef=='None':
            # no temperature cost to be considered
            self.net_dict['T'] = 0
        else:
            raise ValueError('T dependence of '+self.Tdef+' not recognised!')



    def get_P_pH(self):
        """ Calculate the power cost related to pH and update
        the net_dict. Make sure you have set pHdef to the defence that
        you want to compute!
        """
        pH_ad_calc =pHadaptations(self.host)
        if self.pHdef=='FluxPerm':
            #use the FluxPerm equation
            self.net_dict['pH']= pH_ad_calc.getFluxPerm_MP()
        elif self.pHdef=='None':
            #no pH cost to be considered
            self.net_dict['pH']=0
        else:
            raise ValueError('pH dependence of '+self.pHdef+' not recognized!')



    def get_P_store(self):
        """Get the net power stored if there is any"""
        if 'store' in self.net_dict:
            return self.net_dict['store']
        else:
            return 0.0


    def compute_P_growth(self, P_supply):
        """Compute and return the power that can go into growing new biomass
        """

        self.update_P_loss(P_supply)
        return (P_supply * (1.0 - self.P_loss))


    def get_netdictstr(self):
        """Get the net maintenance dictionary as a sring"""
        ###################################

        netdictstr='{'
        for key in sorted(self.net_dict):
            netdictstr += "'" +key + "' : " + str(self.net_dict[key]) +' , '
            # compdict[key] = E.composition[key].activity
        netdictstr+='}'
        return netdictstr

    def set_from_netdictstr(self, netdictstr):
        """Set the net maintenance dicitonary from a string"""
        self.net_dict = ast.literal_eval(netdictstr)
