
import numpy as np
import csv
from copy import copy, deepcopy

import NutMEG
from NutMEG.util.loggersetup import loggersetup as logset
import os
import time as timer
from .ecosystem_dbhelper import db_helper
import math

import NutMEG.util.NutMEGparams as nmp
from NutMEG.util.loggersetup import loggersetup as logset
logger = logset.get_logger(__name__, filelevel=nmp.filelevel, printlevel=nmp.printlevel)

class ecosystem:
    """
    Module to manage a unique ecosystem, containing a ``culture``
    of organisms and their local environment (``reactor``)

    Attributes
    ----------
    r : ``reactor`` like
        Environment information and calculations
    c : ``culture``
        Organism information and calculations
    backup_r : ``reactor`` like
        duplicate of the initial ``r``
    stoppingdict : dict
        Dictionary to tell ecosystem when to stop. See ``setup_stoppingdict``
        for details
    SimID : str
        Simluation identifier for use in data stroage
    dbh : ``ecosystem_dbhelper``
        Helper attribute for datbase management.
    """

    r = None # reactor
    c = None # culture


    def __init__(self, r, c, loggerlevel='INFO', dbpath=nmp.std_dbpath):
        self.r = r
        self.c = c
        self.backup_r = deepcopy(r)
        self.dbh = db_helper(self, dbpath)
        self.stoppingdict=self.setup_stoppingdict()

    def setup_stoppingdict(self):
        """ Initialise the default stoppingdict to be used to tell ecosystem
        when to stop. Each entry has another dict as a value with keys 'Max',
        'Min, 'Consistency', 'Count'. 'Max' and 'Min' are the maximum and
        minimum values of some parameter, 'Consistency' is ho many time steps
        in a row the ecosystem must be outside this limits to stop, and
        'Counter' is the a rolling count of that number.

        Returns
        -------
        The default stoppingdict. Can be edited by updating
        ``ecosystem.stoppingdict``
        """
        stVolume_Fraction = {'Max':0.99, 'Min':0., 'Consistency':0, 'Count':0}
        stMaintenance_Fraction = {'Max':1.0, 'Min':-0.1, 'Consistency':10, 'Count':0}
        stMetabolic_Rate = {'Max':float('inf'), 'Min':1e-40, 'Consistency':10, 'Count':0}
        stGrowth_Rate = {'Max':float('inf'), 'Min':-0.5, 'Consistency':50, 'Count':0}
        stPopulation = {'Max':float('inf'), 'Min':0, 'Consistency':0, 'Count':0}

        return {'Volume_Fraction' : stVolume_Fraction,
          'Maintenance_Fraction' : stMaintenance_Fraction,
          'Metabolic_Rate': stMetabolic_Rate,
          'Growth_Rate' : stGrowth_Rate,
          'Population' : stPopulation}


    def predict_growth(self, tmax=1e17, dt=None, factorup=1.01,
      staticdt=False, quiet=False):
        """Predict how a single culture will evolve in a single environment.

        Simulate the culture growing at a self-determined timestep (unless
        otherwise specified up until
        either a limiter of ``stoppingdict`` is reached, ``tmax`` is exceeded,
        or all of the organisms become inactive.

        Parameters
        ----------
        tmax : float, optional
            Maximum time in s to simulate. Default is 1e17 (3 billion years...)
        dt : float, optional
            If you want to specify the timestep, pass it as dt in s. If dt is
            not passed, one will be selected automatically
        factorup : float, optional
            The factor with which you'd like the culture to grow by when
            automaticall selecting the timestep. Default 1.01
        staticdt : bool, optional
            Pass as ``True`` if you don't want the timestep to be dynamically
            updated every 10000 steps (may result in code running very slow
            for long simulations). Default False.
        quiet : bool, optional
            By default current parameters are printed out to the terminal every
            100 steps. To stop this, pass quiet as ``True``. Default False.
        """

        resultsdict = self.dbh.buildresultsdict()
        logger.info('Checking if this is a new simulation')
        is_new_sim, self.SimID = self.dbh.createtable()

        if not is_new_sim:
            abortstr = ('This combination of organism(s) and location has ' + \
            'been performed before. \n \n Aborting simulation, for results ' + \
            'find: \n \t SimID: '+self.SimID)
            logger.info(abortstr)
            if not quiet:
                print(abortstr)
            return


        logger.info('Setting up growth prediction.')
        stepmax = 100
        t=0
        if dt==None:
            dt = self.c.getmin_timestep(factorup=factorup)
            logger.debug('\n Simulation timestep selected: '+ str(dt))
        else:
            staticdt=True

        self.stoppingdict['Growth_Rate']['Min'] = (
          self.stoppingdict['Growth_Rate']['Min']/dt)

        first=True # this is the first step
        timestart=timer.time() # monitor how long simulations take in real time
        step=1
        while t < tmax:
            timestartstep=timer.time()
            step = step + 1

            logger.debug('\n Starting step '+ str(step))

            startpop = [o.get_population() for o in self.c.all()]
            try:
                self.c.take_step(dt)
            except ZeroDivisionError as zde:
                logger.warning('Zero Division Error: Something may be out of food.')
                logger.info('A Zero Division Error was encountered. Usually' + \
                  'this means that there is no power supply because the '+ \
                  'organism is trying to consume more than is available, ' + \
                  'maintenance fraction is above 1 in the first step, or ' + \
                  "reactions and reagents aren't defined properly")
                if first:
                    logger.info('First step error, saving initial conditions ' + \
                      'as results.')
                    for col in self.c.all():
                        col.output.appendvals(col.get_population(), dt)
                    self.c.output.appendvals()
                resultsdict['Time'].append(t)
                resultsdict['Composition'].append(self.r.getconcs())
                full_results = self.dbh.dict_to_db(resultsdict, end=True)
                break
            except Exception as e:
                logger.info('Other error occurred while taking step')
                logger.info('Saving data up to error...')
                for col in self.c.all():
                    col.output.appendvals(col.get_population(), dt)
                self.c.output.appendvals()
                resultsdict['Time'].append(t)
                resultsdict['Composition'].append(self.r.getconcs())
                full_results = self.dbh.dict_to_db(resultsdict, end=True)
                raise e
            except e:
                raise e

            self.r.update_composition(dt)

            first=False # we made it past the first step!

            t+=dt
            logger.info('Step '+ str(step) + ' complete\n')

            endpop = [o.get_population() for o in self.c.all()]

            resultsdict['Time'].append(t)
            resultsdict['Composition'].append(self.r.getconcs())

            # if sum(endpop) > 0:

            if sum(endpop) <= 0.:
                # full_results = self.dbh.dict_to_db(resultsdict, end=True)
                self.output_step_data(step, full_results)
                logger.info('All of the organisms have died')
                break


            if self.check_stopper('Maintenance_Fraction', self.c.get_maintenance_fractions()):
                break
            if self.check_stopper('Metabolic_Rate', self.c.get_metabolic_rates()):
                break
            if self.check_stopper('Growth_Rate', [abs(g) for g in self.c.get_growth_rates()]):
                break
            if self.check_stopper('Volume_Fraction', [self.c.get_total_volume(inactive=True)]):
                break
            if self.check_stopper('Population', [self.c.get_population(inactive=False)]):
                break


            # but if none of these occur, output step data+save every 100 steps
            if (step/stepmax).is_integer():
                # update on where we are, and save data

                logger.info('Saving Data')
                if t>=tmax:
                    full_results = self.dbh.dict_to_db(resultsdict, end=True)
                else:
                    full_results = self.dbh.dict_to_db(resultsdict, end=False)
                if not quiet:
                    self.output_step_data(step, full_results)

                resultsdict = self.dbh.buildresultsdict() #reset results dictionary

                if (step/10000).is_integer() and not staticdt:
                    # significant local changes are likely, adapt the timestep to compensate
                    dt = self.c.getmin_timestep(factorup=factorup)
                    logger.debug('\n Simulation timestep selected: '+ str(dt))
                    self.stoppingdict['Growth_Rate']['Min'] = (
                      self.stoppingdict['Growth_Rate']['Min']/dt)



        logger.info('Simulation Complete')
        ex_time = timer.time()-timestart
        logger.info('This took '+str(ex_time)+' s  in real time')
        logger.info('Saving Data')
        try:
            full_results = self.dbh.dict_to_db(resultsdict, end=True)
        except:
            logger.warning('Could not save the final state of this simulation!')
        if not quiet:
            self.output_step_data(step, full_results)


    def check_stopper(self, stopper_name, orglst):
        """Check if stoppingdict's requirements are met and we should stop the
        simulation. Returns ``True`` if we should, ``False`` if we shouldn't.

        Parameters
        ----------
        stopper_name : str
            Which stoppper to check (which key in stoppingdict)
        orglst : list
            List of the relevant parameters for each horde/colony in ``c``.

        """
        if min(orglst) > self.stoppingdict[stopper_name]['Max'] or max(orglst) < self.stoppingdict[stopper_name]['Min']:
            self.stoppingdict[stopper_name]['Count'] += 1
            if self.stoppingdict[stopper_name]['Count'] >= self.stoppingdict[stopper_name]['Consistency']:
                logger.info('\n Stopped by ' + stopper_name+'.')
                return True
        else:
            self.stoppingdict[stopper_name]['Count'] = 0
            return False


    def output_step_data(self, step, dicti):
        """Print out key results.

        Parameters
        ----------
        step : int
            Which step of the simulation we are at.
        dicti : dict
            Dictionary of results to show.
        """
        print('\n'+'step = ' + str(step))

        for key in sorted(dicti):
            print('---  '+key+'   '+str(dicti[key][-1]))
