import numpy as np

import NutMEG
# import NutMEG.base_organism as org
import NutMEG.reaction as rxn
from .colony_output import colony_output

import NutMEG.util.NutMEGparams as nmp
from NutMEG.util.loggersetup import loggersetup as logset
logger = logset.get_logger(__name__, filelevel=nmp.filelevel, printlevel=nmp.printlevel)

org = NutMEG.base_organism

class colony():
    """Class for a colony of organisms, which monitors its size and structure.
    """

    collection = np.array([]) # a list of all the active organisms in this colony
    inactive = np.array([]) # a list of the inactive organisms
    base = None # the base organism in the colony, used to make new ones

    def __init__(self, name, collection=[]):
        self.name = name
        self.collection=np.array(collection)
        if len(self.collection) != 0:
            self.base = collection[0].reproduce()
            self.OrgID=self.base.OrgID
            self.output = colony_output(self)



    def add_organism(self, o, n=1):
        """Add an organism or list of organisms to the colony.

        If n is passed, it is the number of random instances of org to
        generate.
        """
        # if this is the first organism to join the colony, save it as
        # the base.
        if len(self.collection) == 0:
            self.base = o.reproduce()
            self.OrgID=self.base.OrgID
            self.output = colony_output(self)

        if isinstance(o, org):# is org.organism:
            self.collection = np.append(self.collection, o)
        elif type(o) is list:
            for a in o:
                if isinstance(o, org):
                        self.collection = np.append(self.collection, o)
                else:
                    raise TypeError(str(o) + 'is not an organism!')
        else:
            raise TypeError('You have not passed an organism or list of ' + \
              'organisms!')


        if n!=1:
            new_population=[]
            for i in range(1, n):
                # prepare the base organism for splitting
                self.base.E_growth += self.base.E_synth
                self.base.volume += self.base.base_volume
                new = self.base.reproduce()

                # let them have random levels of growth
                new.E_growth = (new.E_synth*np.random.rand())
                new_population.append(new)
            self.collection = np.append(self.collection, new_population)


    def take_step(self, t):
        """ Advance the colony by time t.

        We calculate the local environment only once to save doing it
        for every organsim.
        """

        self.base_Ps = self.base.get_supplied_power(update_energetics=True)

        deceased = [] # List the organisms that must be removed
        counter = 0
        startnum = len(self.collection)
        for o in self.collection:
            if o.isactive:
                o.respiration.net_pathway.molar_gibbs = ( \
                  self.base.respiration.net_pathway.molar_gibbs)
                o.respiration.G_C = self.base.respiration.G_C
                o.metabolic_rate = self.base.metabolic_rate
                o.respiration.rate = self.base.metabolic_rate
                o.take_step(t, update_energetics=False)
                if o.issplitting:
                    # Add in a freshly spawned, independent organism.
                    self.add_organism(o.reproduce())
            else:
                # the organism is dormant, move it out of the active collection
                self.inactive = np.append(self.inactive, o)
                deceased.append(counter)
            counter += 1
        # remove the now inactive organisms
        np.delete(self.collection, deceased)


        # update the colony with its new composition?
        # this could be doing nothing, and could also accidentally be replenishing stuff. Check, and debug if necessary.
        self.base.locale.unify_reaction(
          self.base.locale.reactionlist[
            self.base.respiration.net_pathway.equation][
              type(self.base.respiration.net_pathway)],
          overwrite=True)

        #update the data output
        self.output.appendvals(startnum, t)


    def select_timestep(self, factorup=1.01):

        dt = 0.005/1.2
        Egrow=0.
        while Egrow < (factorup-1)*org.E_synth:

            dt = dt*1.2 # amke the timestep 20% bigger
            logger.debug('Trying: ' + str(dt) + ' ... Prev. Energy retrieved = '
              + str(Eretrieved) + ' J')

            if dt > 365*24*3600*1e9:
                # There is no growth in a billion years
                # pass this as the maximum time step
                return 365*24*3600*1e9

            org = deepcopy(self.collection[0]) # as not to meddle with the org

            org.E_growth=0. # in the colony initialisation they're randomised
              # so knock it back down to 0.


            try:
                org.take_step(dt, update_energetics=True)
                Egrow = org.E_growth
            except:
                logger.debug('Error encountered while trying this timestep')



        logger.info('Min timestep for ' + self.name + ': ' + str(dt) + ' s')

        return dt



    def get_population(self, inactive=False):
        """Return the number of organisms in the colony
        """
        if inactive:
            return(len(self.collection)+len(self.inactive))
        else:
            return(len(self.collection))

    def get_mass(self, inactive=False):
        """Return the total biomass (dead or alive) of the colony.
        """
        biomass = 0.
        if inactive:
            for o in np.append(self.collection, self.inactive):
                biomass += o.mass
            return biomass
        else:
            for o in self.collection:
                biomass += o.mass
            return biomass

    def get_volume(self, inactive=False):
        """Return the total volume (dead or alive) of the colony.
        """
        vol = 0.
        if inactive:
            for o in np.append(self.collection, self.inactive):
                vol += o.volume
            return vol
        else:
            for o in self.collection:
                vol += o.volume
            return vol

    def get_average_maintenance_fraction(self):
        """Return the average maintenance fraction of energetic input
        of the alive population.
        """
        avg = 0.
        for o in self.collection:
            #avg += sum(o.maintenance.frac_dict.values())
            avg += sum(o.maintenance.frac_dict.values())
        avg = avg / self.get_population()
        return avg

    def getCHNOPSut(self):
        """ Return the rate of uptake of each CHNOPS element for the
        whole colony.

        When CHNOPSexchanger gets an overhaul, this will need to be
        adapted.
        """
        utlst = []
        for col in self.collection:
            utlst.append(col.CHNOPS.get_uptake())
        return utlst
