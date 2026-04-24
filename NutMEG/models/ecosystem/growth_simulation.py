import numpy as np
import ast
from copy import deepcopy
from NutMEG.core.resident.organism_population import OrganismPopulation

class GrowthSimulation():

    def __init__(self, ES, stoppingdict='default', resultsdict=None):
        self.ES = ES
        if stoppingdict == 'default':
            self.setup_stoppingdict()
        else:
            self.stoppingdict = stoppingdict

        self.resultsdict = {}
        for o in self.ES.residents:
            if isinstance(o, OrganismPopulation):
                self.resultsdict[o.base.name+'_num'] = o.get_num

        if resultsdict:
            # user needs to pass a dict of the form {key : func}, where func
            # will retrieve the attribute you want.
            self.resultsdict.update(resultsdict)


    def run(self, dt, tmax):

        # set up resultsdict
        num_steps = int(tmax / dt)
        _thisrun = {'t':np.zeros(num_steps)}
        for k,v in self.resultsdict.items():
            _thisrun[k] = np.zeros(num_steps)

        for i in range(num_steps):

            self.ES.take_step(dt)

            _thisrun['t'][i] = i*dt
            for k,v in self.resultsdict.items():
                # extract requested results
                _thisrun[k][i] = float(v())

        return _thisrun


    def setup_stoppingdict(self):
        """ Initialise the default stoppingdict to be used to tell the simulation
        when to stop. Each entry has another dict as a value with keys 'Max',
        'Min, 'Consistency', 'Count'. 'Max' and 'Min' are the maximum and
        minimum values of some parameter, 'Consistency' is how many time steps
        in a row the ecosystem must be outside this limits to stop, and
        'Counter' is a rolling count of that number.

        Returns
        -------
        The default stoppingdict. Can be edited by updating
        ``GrowthSimulation.stoppingdict``
        """
        stVolume_Fraction = {'Max':0.99, 'Min':0., 'Consistency':0, 'Count':0}
        stMaintenance_Fraction = {'Max':1.0, 'Min':-0.1, 'Consistency':10, 'Count':0}
        stMetabolic_Rate = {'Max':float('inf'), 'Min':1e-40, 'Consistency':10, 'Count':0}
        stGrowth_Rate = {'Max':float('inf'), 'Min':-0.5, 'Consistency':50, 'Count':0}
        stPopulation = {'Max':float('inf'), 'Min':0, 'Consistency':10, 'Count':0}
