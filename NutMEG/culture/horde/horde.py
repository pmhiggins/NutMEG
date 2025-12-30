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
      unit='cells',
      deathrate=0.,
      biomass_cell_ratio=1.5,
      workoutID=True,
      **bo_kwargs):

        self.name = name
        self.num = num
        self.unit = unit
        self.deathnum=0.

        # don't let the parent initialisation set the ID just yet
        wID = bo_kwargs.pop('workoutID', workoutID)
        bo_kwargs['workoutID'] = False

        super().__init__(name, locale, metabolism, **bo_kwargs)

        self.deathrate = deathrate
        if self.base_life_span and self.deathrate != 0:
            raise ValueError('Horde initialised with a base_life_span and ' +\
              'a deathrate! Please choose one or the other!')
        if self.deathrate < 0:
            raise ValueError('Negative death rate is unphysical')

        # self.maintenance = maintainer(self,
        #   Tdef=kwargs.pop('Tdef', 'None'), pHdef=kwargs.pop('pHdef', 'None'),
        #   Basal=kwargs.pop('Basal',0.0))

        self.biomass_cell_ratio=biomass_cell_ratio
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

        self.growth_rate = (change / (self.num*t))

        self.num += change
        self.volume += new_biomass_cells*self.base_volume
        if change <0:
            self.deathnum += -change


        if self.base_life_span:
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
                self.volume=0.





    def take_step(self, t, update_energetics=False):
        """Overwrite base_organisms take_step. Send the horde forward by time t.
        Perform all metabolic reactions and grow the horde if possible.
        """
        logger.debug(self.OrgID + ' taking step.')

        startnum = deepcopy(self.num)
        if startnum == 0:
            self.output.appendvals(t)
            return

        self.age += t

        self.P_s = self.get_supplied_power(update_energetics=update_energetics)
        logger.debug(self.OrgID+' supplied power = ' + str(self.P_s))

        self.P_EL_growth = self.maintenance.compute_P_growth(
          self.P_s)
        logger.debug(self.OrgID+' EL growth power = ' + str(self.P_EL_growth))

        self.E_store += self.maintenance.get_P_store()*self.num*t
        logger.debug(self.OrgID+' energy store = ' + str(self.E_store))

        P_G_net = self.CHNOPS.grow_with_nutrients(t, numcells=self.num)
        # net power used for growth is passed back.
        # there is a chance this can be negative, if denaturation and
        # nutrient limitation are important.

        # if P_G_net is greater than 0, it corresponds to all the power
        # from P_growth that can acually go into growing new biomass.
        # If there is a maintenance process that requires rebuilding, that
        # energy cost was already accounted for in the respirator.
        if P_G_net > 0.:
            self.P_growth = P_G_net # energy and nutrient limited growth power
            self.P_s -= (self.P_EL_growth - self.P_growth)

        # if P_G_net is less than 0, nutrient availability prevents biomass
        # synthesis so strongly that no energy will be useful for new growth.
        # As the repair energy cost is already factored in, the maximum
        # useful P_S is equal to P_M and no (or negative) growth occurs.
        elif P_G_net <= 0.:
            self.P_growth =P_G_net
            self.P_s -= self.P_EL_growth


        logger.debug(self.OrgID+' net growth power = ' + str(self.P_growth))

        # reduce the cell specific respiration rate accordingly
        self.respiration.rate = self.P_s / self.respiration.G_C

        # perform the catabolic reaction with the locale
        moles_consumed = self.num*self.respiration.rate*t
        self.locale.perform_reaction(self.respiration.net_pathway.equation,
          moles_consumed, re_type=type(self.respiration.net_pathway))

        self.E_growth = self.P_growth * t # instantaneous cell-specific energy into growth

        E_growth_step = self.P_growth*t*self.num # the amount of energy
          # going into growth this step for the whole horde

        new_biomass_cells = E_growth_step/self.E_synth
        self.update_num_vol(t, new_biomass_cells)

        self.output.appendvals(t)
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
