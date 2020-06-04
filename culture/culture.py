"""
Module to manage all organism behaviour in NutMEG. A culture can contain
one or more hordes and colonies, each of those based on (or including) a
base_organism. Each horde/organism requires a maintainer, respirator and
CHNOPSexchanger object (or, they can just use the default).
"""


import NutMEG
import numpy as np
from itertools import chain
from .culture_output import culture_output

import NutMEG.util.NutMEGparams as nmp
from NutMEG.util.loggersetup import loggersetup as logset
logger = logset.get_logger(__name__, filelevel=nmp.filelevel, printlevel=nmp.printlevel)

class culture:
    """Class to contain all organisms in the ecosystem, separated into
    hordes and colonies.

    Attributes
    ----------
    hordes : lst
        List of different hordes in the system
    colonies : lst
        List of different colonies in the system
    lite : bool, optional
        Whether the colonies are `lite` or not. lite colonies are not yet
        supported so please leave this as the default (False)
    output : culture_output
        Helper instance for generating the output and database filling info.
    """

    lite = False # set True for colony_lites, should we include one.


    def __init__(self, hordes=[], colonies=[], lite=False):
        self.hordes = np.array(hordes)
        self.colonies = np.array(colonies)
        self.lite = lite

        self.output = culture_output(self)


    def take_step(self, t):
        """Advance all organisms by the passed time step"""
        for h in self.hordes:
            h.take_step(t)
        if not self.lite:
            for c in self.colonies:
                c.take_step(t)
        # add in else should we include colony_lite.

        self.output.appendvals()

    def getmin_timestep(self, factorup=1.01):
        """Return the smallest time step for all of the organisms."""
        dtlst = []
        for h in chain(self.hordes, self.colonies):
            dtlst.append(h.select_timestep(factorup=factorup))
        return min(dtlst)


    # PARAMETER EXTRACTION

    def get_population(self, inactive=False):
        """Return the total number of hosted organisms"""
        num = 0.
        for o in self.all():
            num += o.get_population(inactive=inactive)
        return num

    def get_total_mass(self, inactive=False):
        """Return the total mass of the hosted organisms"""
        TM = 0.
        for o in self.all():
            TM += o.get_mass(inactive=inactive)
        return TM

    def get_total_volume(self, inactive=False):
        """Return the total volume of the hosted organisms"""
        TV = 0.
        for o in self.all():
            TV += o.get_volume(inactive=inactive)
        return TV



    def get_metabolic_rates(self):
        """ Return a list of the metabolic rates of the hosted organisms in the
        order hordes, colonies."""
        mr = []
        for h in self.hordes:
            mr.append(h.respiration.rate)
        if not self.lite:
            for c in self.colonies:
                mr.append(c.collection[0].respiration.rate)
        return mr

    def get_growth_rates(self):
        """ Return a list of the growth rates of the hosted organisms in the
        order hordes, colonies."""
        gr = []
        for h in self.hordes:
            gr.append(h.output.params['GrowthRate_'+h.OrgID][-1])
        if not self.lite:
            for c in self.colonies:
                gr.append(c.output.params['GrowthRate_'+c.OrgID][-1])
        return gr

    def get_maintenance_fractions(self):
        """ Return a list of the maintenance fractions of the hosted organisms
        in the order hordes, colonies."""
        mf = []
        for h in self.hordes:
            mf.append(h.output.params['MaintenanceFrac_'+h.OrgID][-1])
        if not self.lite:
            for c in self.colonies:
                mf.append(c.output.params['MaintenanceFrac_'+c.OrgID][-1])
        return mf


    def all(self):
        """Return a list of hosted organisms in the order, hordes, colonies."""
        return chain(self.hordes, self.colonies)

    def full_output(self):
        """generate a full dictionary of everything being monitored in the
        culture, for every colony and horde.
        """
        out = {}
        out.update(self.output.params)
        for h in self.hordes:
            out.update(h.output.params)
        if not self.lite:
            for c in self.colonies:
                out.update(c.output.params)

        # add in else should we include colony_lite.
        # print(out['totBM_cells'])
        # print(out['no_alive_SlowHorde'])
        return out

    def refresh_output(self):
        """Refresh all of the data output dictionaries in the culture"""
        for h in self.all():
            h.output.refreshparams()
        self.output.refreshparams()


    def full_output_datatypes(self):
        """Return the name and type of all parameters output monitors. Useful
        for database building."""
        out = self.output.paramtypelst()
        for h in self.hordes:
            for op in h.output.paramtypelst():
                out.append(op)
        if not self.lite:
            for c in self.colonies:
                for op in c.output.paramtypelst():
                    out.append(op)
        return out
