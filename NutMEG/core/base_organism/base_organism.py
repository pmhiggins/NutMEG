# import sys, os, math, ast
# sys.path.append("../..")
# from copy import deepcopy
# import NutMEG.reaction as rxn
# from NutMEG.environment import environment
# from NutMEG.core.reactor import reactor
# from .maintainer import maintainer
# from .CHNOPSexchanger import CHNOPSexchanger
# from .respirator import respirator
# from .base_organism_dbhelper import bodb_helper
# from .synthesis.cell_synthesis import cell_synthesis as synth
import math
from .grower import Grower
from .maintainer import Maintainer
from .metaboliser import Metaboliser
from NutMEG.core.reactor.reaction import Reaction


# import NutMEG.util.NutMEGparams as nmp
# from NutMEG.util.loggersetup import loggersetup as logset
# logger = logset.get_logger(__name__, filelevel=nmp.filelevel, printlevel=nmp.printlevel)

def SAget(v):
    return (4*math.pi*((v*3)**2))**(1/3)

class BaseOrganism:
    """
    This is a parent class for a model organism, to be used as the base
    reference point for organismic behaviour. It has a large number of attributes,
    but many are optional. The default is usually to behave like E. Coli.

    Attributes
    ----------
    name : str
        name of model organism, used for ecosystem management and output.
    metabolism : metaboliser
        manage all metabolic behaviour such as rate of substrate uptake.
    maintenance : maintainer, optional
        manage all maintnenace behaviour (Default is None, if None is passed on
        initiation the organism will have a maintainer with no survial costs)
    growth : grower, optional
        Manages growth rates. Default behaviour is to use a Bioenergetic
        growth model.
    mass : float, optional
        Mass of each individual cell in kg. Default 1e-15
    dry_mass : float
        Dry mass of each individual cell in kg (e.g. its mass after being
        dehydrated) this is typically around 30% the mass. Default 3e-16
    volume : float
        volume of each individual cell in m^3. Default 1e-18
    age : float, kwarg
        Age of the organism in seconds. Default 0.
    E_synth : float, kwarg
        The amount of energy required to synthesise one cell in the locale.
        Default is 8x10^(-10) J (Higgins 2022). In NutMEG v1, it is possible
        to estimate using the synthesis module. For v2+ this is TBC.
    surfacearea : float, kwarg
        surface area of the cell im m^2 Default is to caluclate from volume
        assuming a spherical organism.
    base_life_span : float, kwarg
        Maximum age for the organism before it becomes inactive in s. Default
        is infinite.
    protein_fraction : float, kwarg
        Fraction of the organisms dry mass which is made up of proteins.
        Useful for synthesis calculations. Default is 0.55 (E. Coli)
    isactive : bool
        Whether the organism is active, e.g. interacting with the locale,
        metabolising, growing, etc. Starts as True, Changed to False when
        neccessary.
    """


    def __init__(self, name, metabolism,
      growth=None,
      maintenance=None,
      mass=1.e-15,
      dry_mass=3e-16,
      volume=1.e-18,
      surfacearea=None,
      E_synth=8e-10,#was 8e-10
      age=0.,
      base_life_span=None):
        """
        Parameters
        ----------
        name : str
            name of model organism, used for ecosystem management and output.
        metabolism : metaboliser
            manage all metabolic behaviour such as rate of substrate uptake.
        maintenance : maintainer, optional
            manage all maintnenace behaviour (Default is None, if None is passed on
            initiation the organism will have a maintainer with no survial costs)
        growth : grower, optional
            Manages growth rates. Default behaviour is to use a Bioenergetic
            growth model.
        mass : float, optional
            Mass of each individual cell in kg. Default 1e-15
        dry_mass : float
            Dry mass of each individual cell in kg (e.g. its mass after being
            dehydrated) this is typically around 30% the mass. Default 3e-16
        volume : float
            volume of each individual cell in m^3. Default 1e-18
        age : float, kwarg
            Age of the organism in seconds. Default 0.
        E_synth : float, kwarg
            The amount of energy required to synthesise one cell in the locale.
            Default is 8x10^(-10) J (Higgins 2022). In NutMEG v1, it is possible
            to estimate using the synthesis module. For v2+ this is TBC.
        surfacearea : float, kwarg
            surface area of the cell im m^2 Default is to caluclate from volume
            assuming a spherical organism.
        base_life_span : float, kwarg
            Maximum age for the organism before it becomes inactive in s. Default
            is infinite.
        protein_fraction : float, kwarg
            Fraction of the organisms dry mass which is made up of proteins.
            Useful for synthesis calculations. Default is 0.55 (E. Coli)
        isactive : bool
            Whether the organism is active, e.g. interacting with the locale,
            metabolising, growing, etc. Starts as True, Changed to False when
            neccessary.
        """

        self.name=name
        self.mass = mass
        self.dry_mass=dry_mass
        self.volume = volume
        self.surfacearea = surfacearea
        if not self.surfacearea:
            self.surfacearea = SAget(self.volume)

        if isinstance(metabolism, Metaboliser):
            self.metabolism = metabolism
        elif isinstance(metabolism, Reaction):
            self.metabolism = Metaboliser(metabolism)
        else:
            raise ValueError("base_organism cannot be initialised without a net metabolic reaction")

        if not maintenance:
            self.maintenance = Maintainer([])
        elif type(maintenance) is Maintainer:
            self.maintenance = maintenance
        else:
            warnings.warn('base_organism: '+name+' has been initalised without a maintainer')

        if not growth:
            self.growth = Grower()
        elif type(growth) is Grower:
            self.growth = growth
        else:
            warnings.warn('base_organism: '+name+' has been initalised without a grower')





        # self.respiration = respirator(self, metabolism, **respiration_kwargs)
        # self.age = age
        # self.mass=mass
        # self.dry_mass=dry_mass
        # self.volume = volume
        # self.base_volume = volume
        self.E_synth = E_synth #was 8e-10

        # self.max_metabolic_rate = kwargs.get('max_metabolic_rate', float('inf'))
        # self.update_metabolic_rate()
        # self.E_store = E_store
        # self.pH_interior = pH_interior
        # self.membrane_potential = membrane_potential
        # self.membrane_permH = membrane_permH
        # self.membrane_permOH = membrane_permOH

        # self.surfacearea = surfacearea
        # if not self.surfacearea:
        #     self.surfacearea = SAget(self.volume)

        # self.E_growth=0.0

        self.isactive=True
        self.issplitting=False
        self.base_life_span = base_life_span

        # if maintenance == None:
        #     self.maintenance = maintainer(self, **maintenance_kwargs)
              # Tdef=kwargs.get('Tdef', 'None'), pHdef=kwargs.get('pHdef', 'None'),
              # net_dict={'Basal':kwargs.get('Basal',0.0)},
              # rebuild=kwargs.get('rebuild',[]))
        # else:
        #     self.maintenance = maintenance
        # if CHNOPS == None:
        #     self.CHNOPS = CHNOPSexchanger(self, **CHNOPS_kwargs)
        # else:
        #     self.CHNOPS = CHNOPS
        #
        # self.dbh = bodb_helper(self, dbpath=dbpath)
        # if workoutID:
        #     self.dbh.workoutID()
        # else:
        #     self.OrgID = self.name
        #
        # self.P_s = self.get_supplied_power()
        # self.P_EL_G = 0.
        # self.P_G = 0.

    def update(self, locale):
        """
        Update the following time- and location-dependent organismic parameters:
        - metabolic rate
        - maintenance requirement (if present)
        - growth rate.
        This function does not update the reactor composition or free energies.
        """
        self.metabolism.compute_rate(self, locale)
        self.maintenance.compute_maintenance(self, locale)
        self.growth.compute_rate(self, locale)


    # commented out functions below need to be updated with the new
    # attributes for base_organism. But do not need to be reviewed until we
    # review the whole database system.

    # @classmethod
    # def bo_from_db(cls, name, locale, OrgID, num=1, dbpath=nmp.std_dbpath):
    #     """Create an instance of a known organism from an entry in a NutMEG
    #     database.
    #
    #     Parameters
    #     ----------
    #     name : str
    #         Name of organism to extract, determines table to use
    #     locale : reactor like
    #         Reactor object for the organism to exist in
    #     OrgID : str
    #         OrgID of the database entry to base this organism from
    #     num : int, optional
    #         number of organisms, for if this is a horde instance
    #     dbpath : str, optional
    #         path of the dictionary to extract from.
    #     """
    #
    #     dbdict = bodb_helper.from_db(name, OrgID, dbpath=dbpath)
    #
    #     #TODO  add in CHNOPS!
    #     O = cls(name, locale, dbdict['Respiration'][1], num=num,
    #       maintenance = None,
    #       workoutID = False,
    #       mass = dbdict['Mass'][1],
    #       dry_mass = dbdict['DryMass'][1],
    #       E_synth = dbdict['Esynth'][1],
    #       volume = dbdict['Volume'][1],
    #       memb_pot = dbdict['MembranePot'][1],
    #       PermH = dbdict['PermH'][1],
    #       PermOH = dbdict['PermOH'][1],
    #       pH_interior = dbdict['pHint'][1],
    #       n_ATP = dbdict['n_ATP'][1],
    #       k_RTP = dbdict['k_RTP'][1],
    #       base_life_span = dbdict['base_life_span'][1])
    #
    #     O.maintenance = maintainer(O,
    #       net_dict=ast.literal_eval(dbdict['MaintenancePower'][1]),
    #       Tdef= dbdict['Tdef'][1], pHdef = dbdict['pHdef'][1])
    #
    #     O.OrgID = OrgID
    #
    #     return O
    #
    #
    # def reset_from_db_dict(self, dbdict):
    #     """Reset this organism's parameters with the ones in dbdict, which must
    #      be a bo_dbhelper output."""
    #     # self.init(self.name, self.locale, self.respirator.net_pathway)
    #     self.E_synth = dbdict['Esynth'][1]
    #     self.dry_mass = dbdict['DryMass'][1]
    #     self.mass = dbdict['Mass'][1]
    #     self.maintenance.set_from_netdictstr(dbdict['MaintenancePower'][1])
    #     self.volume = dbdict['Volume'][1]
    #     self.base_volume = dbdict['Volume'][1]
    #     self.maintenance.Tdef = dbdict['Tdef'][1]
    #     self.maintenance.get_P_T()
    #     self.maintenance.pHdef = dbdict['pHdef'][1]
    #     self.maintenance.get_P_pH()
    #     self.memb_pot = dbdict['MembranePot'][1]
    #     self.PermH = dbdict['PermH'][1]
    #     self.pHinterior = dbdict['pHint'][1]
    #     self.respiration.n_ATP = dbdict['n_ATP'][1]
    #     self.respiration.G_C = self.respiration.n_ATP*self.respiration.G_P
    #     if dbdict['k_RTP'][1] != self.respiration.net_pathway.rate_constant_RTP:
    #         self.respiration.net_pathway.rate_constant_RTP = dbdict['k_RTP'][1]
    #         self.respiration.net_pathway.rate_constant_env = ( \
    #           self.respiration.net_pathway.rate_constant_RTP * \
    #           (2**((self.locale.env.T-298)/10)))


    # @staticmethod
    # def builtin_forcing_funcs(funcID, attrs):
    #
    #     if funcID == 'Monod':
    #         # 2 attrs: substrate ID (e.g., 'H2(aq)'), and Monod half-saturation constant
    #         return (lambda _org, S, K: _org.locale.composition[S].conc/(_org.locale.composition[S].conc + K), [attrs['S'], attrs['K']])
    #     if funcID == 'MineralGoethite':
    #         return (lambda _org, K: (_org.bm_conc/_org.locale.composition['Goethite'].conc)/((_org.bm_conc/_org.locale.composition['Goethite'].conc) + K)), [attrs['K']]
    #
    #     else:
    #         raise ValueError('Unknown custom forcing function bassed to builtin_forcing_funcs')


    # def get_ESynth(self, AA=False, comp=None):
    #     """Use the synthesis module to get the synthesis energy for this
    #     organism. Pass AA as True to include the cost of Amino Acid
    #     synthesis.
    #
    #     By default use E Coli parameters as built into synthesis. A future
    #     update might extend this, meaning we'll have to add comp """
    #
    #     #using comp in this way calls a  database which is already populated
    #     # passing host will only change results byt this organisms' protein_fraction
    #
    #     synth_dict = synth.get_ESynth_density(
    #       self.locale.env.T, AA=AA, compute=comp) #cost of sythesis per dry gram of cells
    #     #update E_synth for this organsim from the calculated density.
    #     if AA:
    #         self.E_synth = synth_dict*self.dry_mass*1000
    #     else:
    #         E_dens = []
    #         [E_dens.append(v) for k,v in synth_dict]
    #         self.E_synth = sum(E_dens)*self.dry_mass*1000



    # def update_metabolic_rate(self):
    #     """
    #     Update the rate of the catablic reaction, without considering
    #     nutrient limitation. Equivalent to calling self.respiration.get_rate()
    #     """
    #     self.respiration.get_rate()
    #     # for the base organism, the metabolic rate is the limiter.
    #     # for other limiters enhance this method to include them
    #     # Perhaps put them in a dictionary/list and find the lowest one.
    #     self.respiration.get_rate()
    #     if self.respiration.rate < self.max_metabolic_rate:
    #         self.metabolic_rate = self.respiration.rate
    #     else:
    #         self.respiration.rate = self.max_metabolic_rate
    #         self.metabolic_rate = self.max_metabolic_rate
    #
    #
    # def get_supplied_power(self, update_energetics=False):
    #     """Using respiration, find the instantaneous power supply by multiplying
    #     the molar gibbs yield and the metabolic rate and return it.
    #
    #     If update_energetics is passed, update the thermodynamic parameters of
    #     the respiration first.
    #     """
    #     if (update_energetics and
    #       (type(self.respiration.net_pathway) is rxn.redox)):
    #         # Use a different pathway, update the reagents first
    #         self.respiration.net_pathway.forward.rto_reagents()
    #         self.respiration.net_pathway.reverse.rto_reagents()
    #         self.respiration.net_pathway.update_E(
    #           getgamma=False, estimateDifferentials=True)
    #         self.respiration.net_pathway.update_molar_gibbs()
    #     elif update_energetics: # just use a normal thermochemical reaction.
    #         # update the energetics of reaction in the current environment
    #         self.respiration.net_pathway.rto_current_env()
    #         # update the free energy of our metabolism
    #         self.respiration.net_pathway.update_molar_gibbs_from_quotient(
    #           updatestdGibbs=False)
    #     self.respiration.get_rate()
    #     logger.debug(self.OrgID+' metabolic rate = ' + str(self.respiration.rate))
    #     if self.respiration.G_C > 0.0 and self.respiration.rate > 0.0:
    #         return (self.respiration.G_C*self.respiration.rate)
    #     else:
    #         # logger.warning(self.OrgID+' has no or negative energy supply!')
    #         return 0.#1e-50


    # def take_step(self, t, update_energetics=False):
    #     """Increment the organism's life by time t.
    #
    #     Updates the organism's parameters based on its mortality,
    #     environment, metabolism, etc. If update_energetics is False, any
    #     thermodynamic parameters in locale will not be updated (ie, T, P dependent ones).
    #     """
    #     if self.isactive:
    #         self.age += t
    #
    #         #get energy flow
    #         self.P_s = self.get_supplied_power(update_energetics)
    #         self.P_EL_growth = self.maintenance.compute_P_growth(self.P_s)
    #
    #         # update energy used in this step
    #         self.E_store += self.maintenance.get_P_store()*t
    #         # E_growth_step = self.P_growth*t # the max amount of energy
    #           # going into growth this step, facilitated by kinetics
    #
    #         P_G_net = self.CHNOPS.grow_with_nutrients(t)
    #           # net power used for growth is passed back.
    #           # there is a chance this can be negative, if denaturation and
    #           # nutrient limitation are important.
    #
    #         # if P_G_net is greater than 0, it corresponds to all the power
    #         # from P_growth that can acually go into growing new biomass.
    #         # If there is a maintenance process that requires rebuilding, that
    #         # energy cost was already accounted for in the respirator.
    #         if P_G_net > 0.:
    #             self.P_growth = P_G_net # energy and nutrient limited growth power
    #             self.P_s -= (self.P_EL_growth - self.P_growth)
    #
    #         # if P_G_net is less than 0, nutrient availability prevents biomass
    #         # synthesis so strongly that no energy will be useful for new growth.
    #         # As the repair energy cost is already factored in, the maximum
    #         # useful P_S is equal to P_M and no (or negative) growth occurs.
    #         else:
    #             # P_G_net <= 0
    #             self.P_growth = 0.
    #             self.P_s -= self.P_EL_growth
    #
    #         # reduce the respiration rate accordingly
    #         self.respiration.rate = self.P_s / self.respiration.G_C
    #
    #         # perform the catabolic reaction with the locale
    #         moles_consumed = self.respiration.rate*t
    #         self.locale.perform_reaction(self.respiration.net_pathway.equation,
    #           moles_consumed, re_type=type(self.respiration.net_pathway))
    #
    #         self.E_growth += self.P_growth * t
    #         # update volume as it grows
    #         self.volume = self.base_volume*(1.0 + (self.E_growth/self.E_synth))
    #
    #         if self.E_growth > self.E_synth:
    #             # split the cell
    #             if round(self.E_growth/self.E_synth)>1:
    #                 raise ValueError('Your timestep is too long, '
    #                   + self.name + ' is splitting multiple times per '
    #                   + 'timestep!')
    #             else:
    #                 self.issplitting = True
    #                 # the organism will be split by the colony.
    #         elif self.E_growth < 0.:
    #             # base_organism shrunk, consider it dead
    #             self.isactive=False
    #
    #
    #     if self.base_life_span:
    #         if self.age > self.base_life_span:
    #             # kill the organism. It remains as biomass,
    #             # but cannot divide or metabolise
    #             self.isactive=False


    # def reproduce(self):
    #     """Return a new organism identical to this one but of age zero.
    #
    #     Also reset the energy stores of the new organism, and ensure its
    #     metabolic reaction is maintained in real time (not a complete
    #     deep copy).
    #     """
    #     new = deepcopy(self)
    #     new.locale = self.locale #make sure everything always points to
    #       #the same locale and not a copy
    #
    #     # to imrove efficiency, give all spawn the parent's respiration object.
    #     new.respiration = self.respiration
    #     new.E_growth=0.
    #     new.E_store=0.
    #     new.age = 0.
    #     new.issplitting = False
    #     self.issplitting = False
    #     self.E_growth -= self.E_synth
    #     self.volume -= self.base_volume
    #     return new


    def __str__(self):
        return self.name
