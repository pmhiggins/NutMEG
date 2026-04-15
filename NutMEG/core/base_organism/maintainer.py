# from .adaptations.pHadaptations import pHadaptations
# from .adaptations.Tadaptations import Tadaptations
# import ast

# import NutMEG.util.NutMEGparams as nmp
# from NutMEG.util.loggersetup import loggersetup as logset
# logger = logset.get_logger(__name__, filelevel=nmp.filelevel, printlevel=nmp.printlevel)

class maintainer:
    """
    This class is for computing and calculating the maintenance
    requirements in the form of powers for a given organism.

    All values are per cell.

    Attributes
    ----------
    mechanisms : list[MaintenanceModel]
        List of MaintenanceModel object representing the different maintenance
        stresses the cell faces.
    total_maintenance_power : float
        Aggregate sum of all maintennce power requirements for host cell.
    MP_list : list
        List same length as ``mechanisms`` detailing the indidual contributions
        to total maintenance of each mechanism modelled.
    rebuild : list, optional
        Which maintenance processes required biomass replacement. List should
        contain string keys corresponding to the maintenance processes that
        require biomass synthesis and hence reduce the energy available for
        new growth. For example, setting ``rebuild = ['T']`` tells the maintainer
        that the energy used to defend against temperature must be used to
        synthesise biomass and hence also consume the organism's nutrient budget.
        Default value is empty (ie, maintenance processes do not need new biomass).


    """

    #P_store = 0. # Fractional power stored away for future use. Depreciated,
    # but I'll reintroduce it if bugs appear.


    def __init__(self, host, locale,
      mechanisms=[]):
      # net_dict={},
      # supply=1.0,
      # Tdef='None',
      # pHdef='None',
      # rebuild=[]):
        """
        Parameters
        ----------
        host : ``base_organism`` like
            host organism.
        locale : ``reactor`` like
            The chemical reactor the organism exists inside.
        mechanisms : list[MaintenanceModel]
            List of MaintenanceModel object representing the different maintenance
            stresses the cell faces.
        """
        self.mechanisms = mechanisms

        self.total_maintenance_power = 0.
        self.total_rebuild_power = 0.

        self.MP_list = [0,]*len(mechanisms)
        # self.host = host # this is a reference to the host organism for
          # e.g. if the environment changes.
        # print('setting '+host.name+' Basal to '+str(Basal))
        # self.net_dict = net_dict
        # self.P_loss = 0.0
        # self.net_dict['Basal']=Basal
        # self.Tdef = Tdef
        # self.get_P_T()
        # self.pHdef = pHdef
        # self.get_P_pH()
        # self.frac_dict={}
        # self.update_frac_dict(supply)
        # self.rebuild = rebuild



    def add_mechanism(self, MM):
        """Add a mechanism to the list of MaintenanceModels"""
        self.mechanisms.append(MM)
        self.MP_list.append(0.)

    def compute_maintenance(self, host, locale):
        _PM = 0.
        for i, MM in enumerate(self.mechanisms):
            this_PM = MM.compute(host,locale)
            self.MP_list[i] = this_PM
            _PM += this_PM
        self.total_maintenance_power = _PM
        return _PM

    # def update_frac_dict(self, P_supply):#, hosthorde=True):
    #     """Update the power loss dictionary by expressing it as a fraction
    #     of the passed power supply.
    #     """
    #     if P_supply != 0:
    #         for key, value in self.net_dict.items():
    #             self.frac_dict[key] = (value/P_supply)
    #     else:
    #         for key, value in self.net_dict.items():
    #             self.frac_dict[key] = (0)
    #         # hacky, stops ValueErrors later in the step
    #         self.host.E_growth=0


    # def update_P_loss(self, P_supply):
    #     """Use the power loss dictionary to get the current total
    #     power lost to maintenance.
    #     """
    #
    #     self.update_frac_dict(P_supply)
    #     loss = 0.
    #
    #     for val in self.frac_dict.values():
    #         loss += (val)#*0.001*random.randrange(500,1500))
    #     if self.host.locale.env.T > (400.): #proteins definitely break down
    #         self.P_loss = 1.0
    #     elif loss<=1.0:
    #         self.P_loss = loss
    #     else:
    #         self.P_loss = 1.0
    #
    #
    #
    # def get_P_T(self):
    #     """ Calculate the power cost related to temperature and update
    #     the net_dict. Make sure you have set Tdef to the defence that
    #     you want to compute! Current options are 'Tijhuis' (Tijhuis et al 1993),
    #     'Lever10pc', 'Lever2pc' (both Lever et al 2015) for 10% and 2%
    #     racemization replacement respectively. Alternatively, 'None' ignores
    #     temperature defenses.
    #     """
    #     T_ad_calc = Tadaptations(self.host)
    #     if self.Tdef=='Lever10pc':
    #         # Use the Lever calculation with replacement at 10% [QE]
    #         self.net_dict['T'] = T_ad_calc.getLeverME()
    #     elif self.Tdef=='Lever2pc':
    #         # Use the Lever calculation with replacement at 2% [QE]
    #         self.net_dict['T'] =  T_ad_calc.getLeverME(cutoff_pc=2)
    #     elif self.Tdef =='Lever1/250':
    #         self.net_dict['T'] =  T_ad_calc.getLeverME(cutoff_pc=0.4)
    #     elif self.Tdef=='Tijhuis':
    #         # use the Tijhuis calculation [QE]
    #         self.net_dict['T'] = T_ad_calc.getTijhuisME()
    #     elif self.Tdef=='TijhuisAerobe':
    #         self.net_dict['T'] = T_ad_calc.getTijhuisAerobe()
    #     elif self.Tdef=='TijhuisAnaerobe':
    #         self.net_dict['T'] = T_ad_calc.getTijhuisAnaerobe()
    #     elif self.Tdef=='TOM':
    #         self.net_dict['T'] = T_ad_calc.getTOM()
    #     elif self.Tdef=='None':
    #         # no temperature cost to be considered
    #         self.net_dict['T'] = 0
    #     else:
    #         raise ValueError('T dependence of '+self.Tdef+' not recognised!')
    #
    #
    #
    # def get_P_pH(self):
    #     """ Calculate the power cost related to pH and update
    #     the net_dict. Make sure you have set pHdef to the defence that
    #     you want to compute! Currently, only 'FluxPerm' is included.
    #     Alternatively 'None' ignores pH defences.
    #     """
    #     pH_ad_calc =pHadaptations(self.host)
    #     if self.pHdef=='FluxPerm':
    #         #use the FluxPerm equation
    #         self.net_dict['pH']= pH_ad_calc.getFluxPerm_MP()
    #     elif self.pHdef=='None':
    #         #no pH cost to be considered
    #         self.net_dict['pH']=0
    #     else:
    #         raise ValueError('pH dependence of '+self.pHdef+' not recognized!')
    #
    #
    #
    # def get_P_store(self):
    #     """Get the net power stored if there is any"""
    #     if 'store' in self.net_dict:
    #         return self.net_dict['store']
    #     else:
    #         return 0.0
    #
    #
    # def compute_P_growth(self, P_supply):
    #     """Compute and return the maximum power that can go into growing new
    #     biomass from the incoming power supply ``P_supply`` in W/cell.
    #     """
    #
    #     self.update_P_loss(P_supply)
    #     return (P_supply * (1.0 - self.P_loss))
    #
    #
    # def calculate_P_rebuild(self):
    #     """
    #     sum together all the maintenance powers that are flagged as
    #     requiring biomass rebuilt
    #     """
    #
    #     self.P_rebuild = 0.
    #     for rb in self.rebuild:
    #         self.P_rebuild += self.net_dict[rb]
    #
    #
    # def get_netdictstr(self):
    #     """Get the net maintenance dictionary as a sring"""
    #     netdictstr='{'
    #     for key in sorted(self.net_dict):
    #         netdictstr += "'" +key + "' : " + str(self.net_dict[key]) +' , '
    #         # compdict[key] = E.composition[key].activity
    #     netdictstr+='}'
    #     return netdictstr
    #
    # def set_from_netdictstr(self, netdictstr):
    #     """Set the net maintenance dicitonary from a string"""
    #     self.net_dict = ast.literal_eval(netdictstr)
