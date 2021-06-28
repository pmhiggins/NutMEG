import NutMEG
from .horde_output import horde_output
from copy import deepcopy
from NutMEG.culture.base_organism.maintainer import maintainer

import NutMEG.util.NutMEGparams as nmp
from NutMEG.util.loggersetup import loggersetup as logset
logger = logset.get_logger(__name__, filelevel=nmp.filelevel, printlevel=nmp.printlevel)


class horde(NutMEG.base_organism):
    """Class for a horde of organisms acting as one for better efficiency.,
    albeit losing some fidelity. For the majority of applications using a horde
    is much better than using a colony. A horde shares most attributes with
    base-organism like objects (and using them as per-cell parameters),
    but has a couple of its own.

    Attributes
    ----------
    num : int
        The number of organisms the horde represents.
    volume : float
        The total volume of the horde.
    deathnum : int
        The number of dead cells which are inactive, but still represent
        biomass.
    biomass_cell_ratio : float, kwarg
        To be used for conversion between cell numbers and volumes, because
        some cells will be mid growth and bigger than others. Default 1.5.
    deathrate : float, kwarg
        Alternative to the lifespan attribute. Implement a death rate which
        represents the fraction of the total horde biomass which becomes
        inactive per second [s^-1]. Default 0.
    """

    def __init__(self, name, locale, metabolism, num,
      maintenance=None,
      CHNOPS=None,
      mass=1e-15,
      dry_mass=3e-16,
      *args, **kwargs):

        self.name = name
        self.num = num
        self.deathnum=0.
        wID = kwargs.pop('workoutID', True)
        kwargs['workoutID'] = False

        NutMEG.base_organism.__init__(self, name, locale,
          metabolism,
          CHNOPS=CHNOPS,
          mass=mass,
          dry_mass=dry_mass,
          *args, **kwargs)

        self.deathrate = kwargs.pop('deathrate', 0.)
        if self.base_life_span < float('inf') and self.deathrate != 0:
            raise ValueError('Horde initialised with a base_life_span and ' +\
              'a deathrate! Please choose one or the other!')
        if self.deathrate < 0:
            raise ValueError('Negative death rate is unphysical')

        # self.maintenance = maintainer(self,
        #   Tdef=kwargs.pop('Tdef', 'None'), pHdef=kwargs.pop('pHdef', 'None'),
        #   Basal=kwargs.pop('Basal',0.0))

        self.biomass_cell_ratio=kwargs.pop('biomass_cell_ratio',1.5)
        self.volume=self.num*self.base_volume*self.biomass_cell_ratio
        self.OrgID = ''
        self.output = horde_output(self)
        if wID:
            self.workoutID()
        self.historicnum=[]
        self.molecons=0 #! counter for how many moles of substrate have been consumed.

    def workoutID(self):
        self.output = horde_output(self)
        self.dbh.workoutID()
        self.output = horde_output(self)

    def reproduce(self):
        """Hordes don't reproduce, throw an error if something tries to
        make it do so.
        """
        raise TypeError('Horde objects cannot reproduce! Did you mean'+\
          'to call update_num_vol?')


    def update_num_vol(self, t, new_biomass_cells):
        """After updating the volume, correct the number of organisms."""

        change = round((self.volume + new_biomass_cells*self.base_volume) / \
          (self.base_volume * self.biomass_cell_ratio)) - self.num

        if self.deathrate > 0.:
            # implement death rate before we add the new biomass
            self.num -= (self.deathrate*self.num*t)
            self.volume -= (self.deathrate*self.volume*t)
            self.deathnum += (self.deathrate*self.num*t)

        self.num += change
        self.volume += new_biomass_cells*self.base_volume

        if self.base_life_span < float('inf'):
            # cells can die, so keep an eye on them.
            self.historicnum.append(change)

            if self.age > self.base_life_span:
                # remove the organisms which have reached their life span.
                lsstep = round((self.age-self.base_life_span)/t) #the step 1 life span ago
                self.num -= self.historicnum[lsstep]
                self.deathnum = self.deathnum + self.historicnum[lsstep]
                self.volume -= (self.historicnum[lsstep] * \
                  (self.base_volume * self.biomass_cell_ratio))
            if self.num <0:
                # make sure we don't go below zero
                self.num=0.
                self.voume=0.





    def take_step(self, t):
        """Overwrite base_organisms take_step. Send the horde forward by time t.
        Perform all metabolic reactions and grow the horde if possible.
        """
        logger.debug(self.OrgID + ' taking step.')

        self.age += t

        self.P_s = self.get_supplied_power(update_energetics=True)
        logger.debug(self.OrgID+' supplied power = ' + str(self.P_s))

        self.P_growth = self.maintenance.compute_P_growth(
          self.P_s)
        logger.debug(self.OrgID+' growth power = ' + str(self.P_growth))

        self.E_store += self.maintenance.get_P_store()*self.num*t
        logger.debug(self.OrgID+' energy store = ' + str(self.E_store))

        E_growth_step = self.P_growth*t*self.num # the amount of energy
          # going into growth this step for the whole horde
        logger.debug(self.OrgID+' net growth energy available = ' + str(E_growth_step))

        if E_growth_step >0:
            #this needs a big check
            E_back = self.CHNOPS.check_nutrients(E_growth_step, t, numcells=self.num)
            # E_back = self.CHNOPS.grow_with_nutrients(E_growth_step, t, numcells=self.num)
        else:
            E_back = 0

        self.E_growth += (E_growth_step - E_back)
        logger.debug(self.OrgID+' total energy to be used for growth = ' + str(self.E_growth))

        startnum = deepcopy(self.num)

        if self.E_growth > 0.:
            self.CHNOPS.grow_with_nutrients(E_growth_step, t, checknutrients=False, ret=E_back)

            new_biomass_cells = self.E_growth/self.E_synth

            # new_biomass_cells, left = divmod(self.E_growth/self.E_synth, 1.0)
            # self.E_growth = left*self.E_synth


            moles_consumed = ((1-(E_back/(self.num*self.P_s*t))) * \
              self.num*self.respiration.rate*t)


            # ((self.E_growth/E_growth_step) * \
            #   self.num*self.respiration.rate*t)
              # The fraction on the front corrects for CHNOPS limitation: the
              # metabolism doesn't need to run so fast.
            logger.debug('Performing '+self.OrgID+"'s reaction with "+str(moles_consumed)+' mol.')
            self.locale.perform_reaction(self.respiration.net_pathway.equation,
              moles_consumed, re_type=type(self.respiration.net_pathway))

            self.E_growth=0
            self.update_num_vol(t, new_biomass_cells)

        else:
            # no growth, but we still need to perform the metabolic reaction.
            moles_consumed = self.respiration.rate*t*self.num

            if self.respiration.net_pathway.molar_gibbs != 0:
                logger.debug('Performing '+self.OrgID+"'s reaction with "+str(moles_consumed)+' mol.')
                self.locale.perform_reaction(self.respiration.net_pathway.equation,
                  moles_consumed, re_type=type(self.respiration.net_pathway))
            else:
                logger.warning('Not enough energy in reactor! Holding '+self.OrgID+'in stasis...')

                self.E_growth = 0
            self.update_num_vol(t, 0.)

        self.output.appendvals(startnum, t)

        self.molecons += moles_consumed



    def select_timestep(self, factorup=1.01, returncop=False):
        """work out a suitable time step for the horde to grow by factorup
        times. return the timestep. If returncop is passed as True, also return
        the copy used to calculate the timestep.
        """
        dt = 0.005/1.2
        cop = deepcopy(self)
        while cop.volume < factorup*self.volume:

            dt = dt*1.2 #make the timestep 20% bigger
            logger.debug('Trying: ' + str(dt) + ' ... Prev. Volume = '
              + str(cop.volume) + ' m^-3')

            if dt > 365*24*3600*1e9:
                # There is no growth in a billion years
                # pass this as the maximum time step
                return 365*24*3600*1e9

            cop = deepcopy(self) # as not to meddle with the organism.

            try:
                cop.take_step(dt)
            except:
                logger.debug('Error encountered while trying this timestep')

        logger.info('Min timestep for ' + self.name + ': ' + str(dt) + ' s')

        if returncop:
            return dt, cop
        else:
            return dt



    def get_mass(self, inactive=False):
        """return the total (approximate) biomass of the horde in kg"""
        if inactive:
            return ((self.num+self.deathnum) * self.mass * \
              self.biomass_cell_ratio)
        else:
            return self.num*self.mass*self.biomass_cell_ratio

    def get_volume(self, inactive=False):
        """return the total volume of the horde in m^3"""
        if inactive:
            #include the volume of inactive biomass
            return self.volume + (self.deathnum * self.base_volume * \
              self.biomass_cell_ratio)
        else:
            return self.volume

    def get_population(self, inactive=False):
        """Return total number of cells in the horde"""
        if inactive:
            return self.num + self.deathnum
        else:
            return self.num
