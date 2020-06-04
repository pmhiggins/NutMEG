import random
import logging

import NutMEG.util.NutMEGparams as nmp
from NutMEG.util.loggersetup import loggersetup as logset
logger = logset.get_logger(__name__, filelevel=nmp.filelevel, printlevel=nmp.printlevel)




class CHNOPSexchanger:
    """Class for managing CHNOPS exchange between an organism and a reactor.

    Attributes
    ----------
    host : ``base_organism`` like
        host organism. Ensure that the host organism's locale object is the
        reactor you want.
    uptake_consts : dict, kwarg
        dictionary of uptake constants for C,H,N,O,P,S etc, in units /s. Default
        is for each to be 10^{-10} (arbitrarily).

    nutrients : dict
        nutrient parameters for each key in C,H,N,O,P,S, there is a list in the
        format [[rate of moles to be used (/s), % needed (for each cell),
        amount in reactor (from reactor.composition),
        max rate constant per cell (from uptake_consts)], [<list of reagents
        which could be used>]
    limter : str
        Name of the most limiting paramters, either a named nutrient or Energy.

    Notes
    -----
    This class is still in a basic (albeit convoluted) state. This is because
    there reamins a lack of literature for this kind of limitation and hopefully
    one day we will be able to inmprove it considerably.
    """

    uptake_consts = {'C':1e-10, 'H':1e-10, 'N':1e-10,
      'O':1e-10, 'P':1e-10, 'S':1e-10}

    def get_default_nutrients(self):
        """set up the nutrients dictionary, using cell composition values of
        E. Coli. This could be adapted in future child classes for specific
        organisms if the data is available.
        """
        # in an ideal world, we would have uptakes and constants from individual
        # chemicals rather than just the elements as shown here.
        return {
          'C':[[0.0, 0.50/12.0, 0.0, self.uptake_consts['C']], []],
          'H':[[0.0, 0.09, 0.0, self.uptake_consts['H']], []],
          'N':[[0.0, 0.15/14.0, 0.0, self.uptake_consts['N']], []],
          'O':[[0.0, 0.20/16.0, 0.0, self.uptake_consts['O']], []],
          'P':[[0.0, 0.04/31.0, 0.0, self.uptake_consts['P']], []],
          'S':[[0.0, 0.02/32.0, 0.0, self.uptake_consts['S']], []]}

      # approx, adapted slightly from E Coli values in microbiology book
      # (Atlas 1995)



    def __init__(self, host, *args, **kwargs):
        self.host = host
        self.limiter = ''
        self.uptake_consts.update(kwargs.pop('uptake_consts', {}))
        self.find_nutrients(init=True)



    def find_nutrients(self, init=False):
        """ Look in the host's locale to find reagents which could work as a
        nutrient source and update nutrients as appropriate.
        Pass init as True to set up the nutrients dict from scratch.

        Notes
        -----
        there must be a way to make this faster. It has to look in the
        reactor every step.
        """

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



    def update_nutrient_yield(self, numcells=1):
        """Get the nutrient concentrations available from the host's locale.

        Parameters
        ----------
        numcells : int, optional
            number of cells we need nutrients for, as this will be greater than
            one for hordes.
        """
        self.find_nutrients()
        MaxMolOrg= {}
        for key, value in self.nutrients.items():
            # work out the max amount of cells we can get from each
            MaxMolOrg[key] = numcells*value[0][2]*value[0][3]/value[0][1]
            # units: dry g / (L s)
        # the smallest value of MaxMolOrg[0] corresponds to the limiting element.
        sublimiter = min(MaxMolOrg, key=MaxMolOrg.get)
        self.limiter=(sublimiter)
        logger.info(self.limiter+' is the limiting nutrient for '+self.host.name)
        # loop over again and assign each with the max we can make
        factor = MaxMolOrg[sublimiter]#/self.nutrients[limiter][0][1]
        for key, value in self.nutrients.items():
            value[0][0] = factor*value[0][1]*(self.host.locale.volume*1000) # convert to L


    def grow_with_nutrients(self, E_growth, t, updatenutrients=True,
      checknutrients=True, ret=0.0, numcells=1):
        """Convert as much of E_growth [J] into biomass as possible in time t
        [s] based on nutrient availability. Remove the nutrients required from
        ``host.locale`` to build that biomass and pass back any leftover energy.


        Parameters
        ----------
        E_growth : float
            Total Energy available for growth in ``t`` in J
        t : float
            Time frame for which this energy is available in s.
        checknutrients : bool, optional
            Whether to check nutreint limitation vs energy limitation before
            analysis. Pass as ``False`` if you want to ignore nutrient
            limitation or set your own ``ret`` for some reason.
            Default ``True``
        updatenutrients : bool, optional
            Whether to update the nutrients dict before checking for nutrients.
            Default True. This only has an effect is checknutrients is also
            ``True``
        ret : float, optional
            The amount of ``E_growth`` to return in J. Only has effect if
            checknutrients is passed as ``False``.
        numcells : int, optional
            Number of cells to be considering for this growth, if a horde is
            being used for example. Default 1.

        Returns
        -------
        E_growth remaining if the host is nutrient limited in J. 0 if the host
        is energy limited.

        Notes
        -----
        This function does NOT grow the ``host`` (e.g. change its volume and
        mass). That should be done elsewhere.
        """
        if checknutrients:
            ret = self.check_nutrients(E_growth, t, updatenutrients=updatenutrients, numcells=numcells)

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
            if food.name != 'H2O(l)':# and food.activity >=1e-12:
                # only remove the fraction we have actually picked up
                # because check_nutrients has throttled E_growth
                self.host.locale.composition[food.name].activity -= (fraction_used*value[0][0]*t)
                self.host.locale.composition[food.name].conc -= (fraction_used*value[0][0]*t)
                #food.activity -= (value[0][0]*t)
        # print(self.host.locale.composition['P(aq)'].activity)
        return ret

    def check_nutrients(self, E_growth, t, updatenutrients=True, numcells=1):
        """Given a timestep ``t`` and potential growth input ``E_growth``,
        compute how much of that energy can be used.

        Parameters
        ----------
        E_growth : float
            Total Energy available for growth in ``t`` in J
        t : float
            Time frame for which this energy is available in s.
        updatenutrients : bool, optional
            Whether to update the nutrients dict before comparing growth
            estimates. Default True.
        numcells : int, optional
            Number of cells to be considering for this growth, if a horde is
            being used for example. Default 1.
        Return
        ------
        Any E_growth which cannot be converted into biomass because of nutrient
        limitation.
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
            # our yields proportionately. This is done in grow_with_nutrients

        elif self.maxg <= 0.0:
            logger.debug('organism(s) fatally substrate limited by '+self.limiter+'!')
            self.host.throttling = 'substrate: '+self.limiter
            return E_growth
        else:
            logger.debug('organism(s) are substrate limited')
            self.host.throttling = 'substrate: '+self.limiter
            # there must be energy left over, send it back
            return E_growth*(1.0-(self.maxg/self.g))



    def print_nutrients(self):
        """Print the nutrient information in a more readable way.
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
