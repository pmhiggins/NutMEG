from .adaptations.pHadaptations import pHadaptations
from .adaptations.Tadaptations import Tadaptations
import ast

import NutMEG.util.NutMEGparams as nmp
from NutMEG.util.loggersetup import loggersetup as logset
logger = logset.get_logger(__name__, filelevel=nmp.filelevel, printlevel=nmp.printlevel)

class maintainer:
    """
    This class is for computing and calculating the maintenance
    requirements in the form of powers for a given organism.

    All values are PER CELL, so for a Horde, when you use these
    attributes you will likely need to scale them up.

    Attributes
    ----------
    host : ``base_organism`` like
        host organism. Ensure that the host organism's locale object is the
        reactor you want.
    net_dict : dict, optional
        dictionary of the net contributions to maintenance in W per cell. By
        default will include contributions of temperature and pH, according to
        Tdef or pHdef. You could also include a 'Basal' maintenance power.
    frac_dict : dict
        dictionary of the contributions to maintenance as a fraction of power
        supply. Shares keys with ``net_dict``.
    Tdef : str, optional
        Which temperature defences to use. Current options are are those from
        Tijhuis (1993), and Lever (2015). See function ``get_P_T`` for options
        which work. For more info, check out the documentation, and if it isn't
        there yet contact me!
    pHdef : str, optional
        Which pH defences to use. See function ``get_P_pH`` for options
        which work. For more info, check out the documentation, and if it isn't
        there yet contact me!
    P_loss : float
        Instantaneous fractional power loss due to all maintenance in W/organism


    """

    #P_store = 0. # Fractional power stored away for future use. Depreciated,
    # but I'll reintroduce it if bugs appear.


    def __init__(self, host,
      net_dict={},
      supply=1.0,
      Tdef='None',
      pHdef='None'):

        self.host = host # this is a reference to the host organism for
          # e.g. if the environment changes.
        # print('setting '+host.name+' Basal to '+str(Basal))
        self.net_dict = net_dict
        self.P_loss = 0.0
        # self.net_dict['Basal']=Basal
        self.Tdef = Tdef
        self.get_P_T()
        self.pHdef = pHdef
        self.get_P_pH()
        self.frac_dict={}
        self.update_frac_dict(supply)



    def add_mechanism(self, name, val):
        """Add a mechanism to the dictionary of coping mechanisms."""
        self.net_dict[name] = val


    def update_frac_dict(self, P_supply):#, hosthorde=True):
        """Update the power loss dictionary by expressing it as a fraction
        of the passed power supply.
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
        you want to compute! Current options are 'Tijhuis' (Tijhuis et al 1993),
        'Lever10pc', 'Lever2pc' (both Lever et al 2015) for 10% and 2%
        racemization replacement respectively. Alternatively, 'None' ignores
        temperature defenses.
        """
        T_ad_calc = Tadaptations(self.host)
        if self.Tdef=='Lever10pc':
            # Use the Lever calculation with replacement at 10% [QE]
            self.net_dict['T'] = T_ad_calc.getLeverME()
        elif self.Tdef=='Lever2pc':
            # Use the Lever calculation with replacement at 2% [QE]
            self.net_dict['T'] =  T_ad_calc.getLeverME(cutoff_pc=2)
        elif self.Tdef =='Lever1/250':
            self.net_dict['T'] =  T_ad_calc.getLeverME(cutoff_pc=0.4)
        elif self.Tdef=='Tijhuis':
            # use the Tijhuis calculation [QE]
            self.net_dict['T'] = T_ad_calc.getTijhuisME()
        elif self.Tdef=='TijhuisAerobe':
            self.net_dict['T'] = T_ad_calc.getTijhuisAerobe()
        elif self.Tdef=='TijhuisAnaerobe':
            self.net_dict['T'] = T_ad_calc.getTijhuisAnaerobe()
        elif self.Tdef=='TOM':
            self.net_dict['T'] = T_ad_calc.getTOM()
        elif self.Tdef=='None':
            # no temperature cost to be considered
            self.net_dict['T'] = 0
        else:
            raise ValueError('T dependence of '+self.Tdef+' not recognised!')



    def get_P_pH(self):
        """ Calculate the power cost related to pH and update
        the net_dict. Make sure you have set pHdef to the defence that
        you want to compute! Currently, only 'FluxPerm' is included.
        Alternatively 'None' ignores pH defences.
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
        from the incoming power supply ``P_supply`` in W/cell.
        """

        self.update_P_loss(P_supply)
        return (P_supply * (1.0 - self.P_loss))


    def get_netdictstr(self):
        """Get the net maintenance dictionary as a sring"""
        netdictstr='{'
        for key in sorted(self.net_dict):
            netdictstr += "'" +key + "' : " + str(self.net_dict[key]) +' , '
            # compdict[key] = E.composition[key].activity
        netdictstr+='}'
        return netdictstr

    def set_from_netdictstr(self, netdictstr):
        """Set the net maintenance dicitonary from a string"""
        self.net_dict = ast.literal_eval(netdictstr)
