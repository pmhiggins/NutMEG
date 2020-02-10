import random
import logging

import NutMEG.util.NutMEGparams as nmp
from NutMEG.util.loggersetup import loggersetup as logset
logger = logset.get_logger(__name__, filelevel=nmp.filelevel, printlevel=nmp.printlevel)

class CHNOPSexchanger:
    """When this gets an overhaul, don't for get to change take_step
    in organism getCHNOPSut in colony
    """

    limiter='' # information as to what the limiter is at the moment

    def get_default_nutrients(self):
        # in an ideal world, we would have uptakes and constants from individual
        # chemicals rather than just the elements as shown here.
        return {
          'C':[[0.0, 0.50/12.0, 0.0, 1e-10], []],
          'H':[[0.0, 0.09, 0.0, 1e-10], []],
          'N':[[0.0, 0.15/14.0, 0.0, 1e-10], []],
          'O':[[0.0, 0.20/16.0, 0.0, 1e-10], []],
          'P':[[0.0, 0.04/31.0, 0.0, 1e-10], []],
          'S':[[0.0, 0.02/32.0, 0.0, 1e-10], []]}
      # in the format [rate of moles to be used, % needed, amount in reactor, max rate constant per cell]
      # followed by a list containing the reactants which can be used
      # approx, adapted slightly from microbiology book



    def __init__(self, host):

        self.host = host
        self.find_nutrients(init=True)


    def find_nutrients(self, init=False):
        # there must be a way to make this faster. It has to look in the
        # reactor every step.
        if init:
            self.nutrients = self.get_default_nutrients()
            for key, value in self.host.locale.composition.items():
                for keyn, valuen in self.nutrients.items():
                    if keyn in key:
                        self.nutrients[keyn][0][2] += value.activity
                        self.nutrients[keyn][1].append(value)
        else:
            for key, value in self.nutrients.items():
                value[0][2] = 0.0
                for r in value[1]:
                    value[0][2] += r.activity

    """
    Not working as expected. I think it has something to do with value[0][1] being the amount of moles per gram. Look at the tester, it seems that when you change the concentration of one (e.g. S, P, the new limiter can have a lower yield than it had before, which doesn't make sense at its concentration hasn't changed. C is the example here.


    EDIT: Think I fixed it but please play some more/fiddle through the maths
    """






    def update_nutrient_yield(self, numcells=1):
        # get the nutrient concentrations available from the reactor
        self.find_nutrients()
        MaxMolOrg= {}
        for key, value in self.nutrients.items():
            # work out the max amount of cells we can get from each
            MaxMolOrg[key] = numcells*value[0][2]*value[0][3]/value[0][1]
            # units: dry g / (L s)
        # the smallest value of MaxMolOrg[0] corresponds to the limiting element.
        sublimiter = min(MaxMolOrg, key=MaxMolOrg.get)
        self.limiter=(sublimiter)
        #print('Looks like '+limiter+' is the limiter')
        # loop over again and assign each with the max we can make
        factor = MaxMolOrg[sublimiter]#/self.nutrients[limiter][0][1]
        for key, value in self.nutrients.items():
            value[0][0] = factor*value[0][1]*(self.host.locale.volume*1000) # convert to L


    def grow_with_nutrients(self, E_growth, t, updatenutrients=True, checknutrients=True, numcells=1):
        """Exchange nutrients with the locale in order to convert as much
        of E_growth into biomass as possible in time t

        Returns what is left of E_growth.
        """
        if checknutrients:
            ret = self.check_nutrients(E_growth, t, updatenutrients=updatenutrients, numcells=numcells)
        else:
            ret = 0.0

        if self.limiter=='Energy':
            # throttle how much we take up to match the incoming growth energy
            for key, value in self.nutrients.items():
                value[0][0] = value[0][0]*self.g/self.maxg

        fraction_used = 1.0-(ret/E_growth)
        for key, value in self.nutrients.items():
            # consume the relevant amount of nutrient
            # where there are multiple sources pick a random one
            comps = [self.host.locale.composition[f.name].activity for f in value[1]]
            food = random.choices(value[1], comps)[0]
            if food.name != 'H2O(l)' and food.activity >=1e-12:
                # only remove the fraction we have actually picked up
                # because check_nutrients has throttled E_growth
                self.host.locale.composition[food.name].activity -= (fraction_used*value[0][0]*t)
                self.host.locale.composition[food.name].conc -= (fraction_used*value[0][0]*t)
                #food.activity -= (value[0][0]*t)
        # print(self.host.locale.composition['P(aq)'].activity)
        return ret

    def check_nutrients(self, E_growth, t, updatenutrients=True, numcells=1):
        """Given a timestep and potential growth input, compute how much of
        that energy can be used.

        Returns the fraction of E_growth which would remain.
        """
        if updatenutrients:
            self.update_nutrient_yield(numcells=numcells)
        # see how many g cell we can possibly make in the time
        # here we use carbon, but any of them should be equivalent it's just a
        # unit change
        self.maxg = t*self.nutrients['C'][0][0]/(self.nutrients['C'][0][1])

        # use E_synth to work out how many grams the energy can yield
        self.g = (E_growth/self.host.E_synth)*self.host.dry_mass*1000 # convert to g

        if self.g < self.maxg:
            logger.debug('organism(s) are energy limited')
            self.limiter='Energy'

            # we can use all the energy in grow_with_nutrients
            return 0.0


            # the energy is limiting, not the substrate, we need to lower
            # our yields proportionately
            # (moved to grow_with_nutrients)
            #for key, value in self.nutrients.items():
            #    value[0][0] = value[0][0]*self.g/self.maxg
        else:
            logger.debug('organism(s) are substrate limited')
            self.host.throttling = 'substrate: '+self.limiter
            # there must be energy left over, send it back
            return E_growth*(1.0-(self.maxg/self.g))



    def print_nutrients(self):
        """Print the nutrient information in a more readable way
        """
        for key, value in self.nutrients.items():
            prettier = []
            for r in value[1]:
                prettier.append(r.name)
            print(key, value[0], prettier)


    def get_uptake(self, numcells=1):
        """ Return a dictionary showing the uptake rate of each nutient"""
        utdictstr='{'
        for key in sorted(self.nutrients):
            utdictstr += ("'" +key+ "': " + \
              str(self.nutrients[key][0][0]/numcells) +', ')
            # compdict[key] = E.composition[key].activity
        utdictstr+='}'
        return utdictstr
