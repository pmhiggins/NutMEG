import sys, os, math, ast
sys.path.append("../..")
from copy import deepcopy
import NutMEG.reaction as rxn
from NutMEG.environment import environment
from NutMEG.reactor import reactor
from .maintainer import maintainer
from .CHNOPSexchanger import CHNOPSexchanger
from .respirator import respirator
from .base_organism_dbhelper import bodb_helper
from .synthesis.cell_synthesis import cell_synthesis as synth


import NutMEG.util.NutMEGparams as nmp
from NutMEG.util.loggersetup import loggersetup as logset
logger = logset.get_logger(__name__, filelevel=nmp.filelevel, printlevel=nmp.printlevel)

def SAget(v):
    return (4*math.pi*((v*3)**2))**(1/3)

class base_organism:
    """ This is a parent class for all organism-like objects which NutMEG
    can use. It's designed to contain everything that both hordes and
    individual organisms have in common.

    """

    name = None
    age = 0
    base_life_span = float('inf')
    issplitting=False
    E_growth = 0
    isactive=True
    proteinfraction=0.55



    def __init__(self, name, locale, metabolism,
      maintenance=None,
      CHNOPS=None,
      workoutID=True,
      mass=1.e-15,
      dry_mass=3e-16,
      volume=1.e-18,
      *args, **kwargs):
        self.name=name
        self.locale = locale
        self.respiration = respirator(self, metabolism,
          kwargs.pop('n_ATP', 1.0), k_RTP=kwargs.pop('k_RTP', None), overwrite=kwargs.pop('overwrite', False))

        self.mass=mass
        self.dry_mass=dry_mass
        self.volume = volume
        self.base_volume = volume
        self.E_synth = kwargs.pop('E_synth', None) #was 8e-10
        if self.E_synth == None:
            self.get_ESynth(AA=True)
        self.update_metabolic_rate()
        self.E_store = kwargs.pop('E_store', 0.)
        self.pH_interior = kwargs.pop('pH_interior', 7.0)
        self.memb_pot = kwargs.pop('memb_pot', 1e-5)
        self.PermH = kwargs.pop('PermH', 1e-10)
        self.PermOH = kwargs.pop('PermOH', 1e-10)
        self.surfacearea = kwargs.pop('surfacearea', SAget(volume))


        self.base_life_span = kwargs.pop('base_life_span', float('inf'))

        if maintenance == None:
            self.maintenance = maintainer(self,
              Tdef=kwargs.pop('Tdef', 'None'), pHdef=kwargs.pop('pHdef', 'None'),
              Basal=kwargs.pop('Basal',0.0))
        else:
            self.maintenance = maintenance
        if CHNOPS == None:
            self.CHNOPS = CHNOPSexchanger(self, *args, **kwargs)
        else:
            self.CHNOPS = CHNOPS

        self.dbh = bodb_helper(self, dbpath=kwargs.pop('dbpath', nmp.std_dbpath))
        if workoutID:
            self.dbh.workoutID()

    @classmethod
    def bo_from_db(cls, name, locale, OrgID, num=1, dbpath=nmp.std_dbpath):

        dbdict = bodb_helper.from_db(name, OrgID, dbpath=dbpath)
        
        # CHNOPS!
        O = cls(name, locale, dbdict['Respiration'][1], num=num,
          maintenance = None,
          workoutID = False,
          mass = dbdict['Mass'][1],
          dry_mass = dbdict['DryMass'][1],
          E_synth = dbdict['Esynth'][1],
          volume = dbdict['Volume'][1],
          memb_pot = dbdict['MembranePot'][1],
          PermH = dbdict['PermH'][1],
          PermOH = dbdict['PermOH'][1],
          pH_interior = dbdict['pHint'][1],
          n_ATP = dbdict['n_ATP'][1],
          k_RTP = dbdict['k_RTP'][1],
          base_life_span = dbdict['base_life_span'][1])

        O.maintenance = maintainer(O,
          net_dict=ast.literal_eval(dbdict['MaintenancePower'][1]),
          Tdef= dbdict['Tdef'][1], pHdef = dbdict['pHdef'][1])

        O.OrgID = OrgID

        return O


    def reset_from_db_dict(self, dbdict):
        """Reset this organism's parameters with ones extracted from
        a database"""
        # self.init(self.name, self.locale, self.respirator.net_pathway)
        self.E_synth = dbdict['Esynth'][1]
        self.dry_mass = dbdict['DryMass'][1]
        self.mass = dbdict['Mass'][1]
        self.maintenance.set_from_netdictstr(dbdict['MaintenancePower'][1])
        self.volume = dbdict['Volume'][1]
        self.base_volume = dbdict['Volume'][1]
        self.maintenance.Tdef = dbdict['Tdef'][1]
        self.maintenance.get_P_T()
        self.maintenance.pHdef = dbdict['pHdef'][1]
        self.maintenance.get_P_pH()
        self.memb_pot = dbdict['MembranePot'][1]
        self.PermH = dbdict['PermH'][1]
        self.pHinterior = dbdict['pHint'][1]
        self.respiration.n_ATP = dbdict['n_ATP'][1]
        self.respiration.G_C = self.respiration.n_ATP*self.respiration.G_P
        if dbdict['k_RTP'][1] != self.respiration.net_pathway.rate_constant_RTP:
            self.respiration.net_pathway.rate_constant_RTP = dbdict['k_RTP'][1]
            self.respiration.net_pathway.rate_constant_env = ( \
              self.respiration.net_pathway.rate_constant_RTP * \
              (2**((self.locale.env.T-298)/10)))


    def get_ESynth(self, AA=True, comp=None):
        """Use the synthesis module to get the synthesis energy for this organism.

        By default use E Coli parameters as built into synthesis. A future update might extend this, meaning we'll have to add comp """

        #using comp in this way calls a  database which is already populated
        # passing host will only change results byt this organisms' proteinfraction

        E_dens = synth.get_ESynth_density(self.locale.env.T, AA=AA, compute=comp) #cost of sythesis per dry gram of cells
        #update E_synth for this organsim from the calculated density.
        self.E_synth = E_dens*self.dry_mass*1000



    def update_metabolic_rate(self):
        """Update the actual rate of material uptake/dismissal, based on the
        slowest processes occuring.
        """
        # for the base organism, the metabolic rate is the limiter.
        # for other limiters enhance this method to include them
        # Perhaps put them in a dictionary/list and find the lowest one.
        self.respiration.get_rate()
        self.metabolic_rate = self.respiration.rate


    def get_supplied_power(self, update_energetics=False):
        if (update_energetics and
          (type(self.respiration.net_pathway) is rxn.redox)):
            # Use a different pathway, update the reagents first
            self.respiration.net_pathway.forward.rto_reagents()
            self.respiration.net_pathway.reverse.rto_reagents()
            self.respiration.net_pathway.update_E(
              getgamma=False, estimateDifferentials=True)
            self.respiration.net_pathway.update_molar_gibbs()
        elif update_energetics: # just use a normal thermochemical reaction.
            # update the energetics of reaction in the current environment
            self.respiration.net_pathway.rto_current_env()
            # update the free energy of our metabolism
            self.respiration.net_pathway.update_molar_gibbs_from_quotient(
              updatestdGibbs=False)
        self.update_metabolic_rate()
        if self.respiration.G_C > 0.0 and self.metabolic_rate > 0.0:
            return (self.respiration.G_C*self.metabolic_rate)
        else:
            logger.warning(self.OrgID+' has no or negative energy supply!')
            return 0.0


    def take_step(self, t, update_energetics=True):
        """Increment the organism's life by time t.

        Updates the organism's parameters based on its mortality,
        environment, metabolism, etc.
        """
        if self.isactive:
            self.age += t

            #get energy flow
            self.P_s = self.get_supplied_power(update_energetics)
            self.P_growth = self.maintenance.compute_P_growth(self.P_s)

            # update energy used in this step
            self.E_store += self.maintenance.get_P_store()*t
            E_growth_step = self.P_growth*t # the max amount of energy
              # going into growth this step, facilitated by kinetics

            E_back = self.CHNOPS.grow_with_nutrients(E_growth_step, t)
              # energy passed back
              # if CHNOPS are limiting growth, less metabolic energy is needed
            self.E_growth += (E_growth_step - E_back)

            # update volume as it grows
            self.volume = self.base_volume*(1.0 + (self.E_growth/self.E_synth))

            # perform the catabolic reaction with the locale
            moles_consumed = (1-(E_back/self.P_s*t) * \
              self.num*self.respiration.rate*t)

            # ((self.E_growth/E_growth_step) * \
            #   self.respiration.rate*t)
            self.locale.perform_reaction(self.respiration.equation,
              moles_consumed, re_type=type(self.respiration))

            if self.E_growth > self.E_synth:
                # split the cell
                if round(self.E_growth/self.E_synth)>1:
                    raise ValueError('Your timestep is too long, '
                      + self.name + ' is splitting multiple times per '
                      + 'timestep!')
                else:
                    self.issplitting = True
                    # the organism will be split by the colony.


        if self.age > self.base_life_span:
            # kill the organism. It remains as biomass,
            # but cannot divide or metabolise
            self.isactive=False


    def reproduce(self):
        """Return a new organism identical to this one but of age zero.

        Also reset the energy stores of the new organism, and ensure its
        metabolic reaction is maintained in real time (not a complete
        deep copy).
        """
        new = deepcopy(self)
        new.locale = self.locale #make sure everything always points to
          #the same locale and not a copy

        # to imrove efficiency, give all spawn the parent's respiration object.
        new.respiration = self.respiration
        new.E_growth=0.
        new.E_store=0.
        new.age = 0.
        new.issplitting = False
        self.issplitting = False
        self.E_growth -= self.E_synth
        self.volume -= self.base_volume
        return new


    def __str__(self):
        return self.name
