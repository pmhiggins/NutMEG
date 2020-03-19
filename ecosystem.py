
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
    Class to manage a unique ecosystem, containing some colony or colonies
    of organisms, and their local environment (reactive system)
    """

    r = None # reactor
    c = None # culture


    def __init__(self, r, c, loggerlevel='INFO', dbpath=nmp.std_dbpath):
        """
        here r is the reactor and c is the culture for brevity.
        """
        self.r = r
        self.c = c
        self.backup_r = deepcopy(r)

        self.dbh = db_helper(self, dbpath)

        self.stoppingdict=self.setup_stoppingdict()


        # logger.__init__(__name__, level=loggerlevel)

    def setup_stoppingdict(self):

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


    def predict_growth(self, tmax=1e17, volmax=0.99, kmin=1e-23,
      max_growthrate=1.0, endMF=True, dt=None, factorup=1.01, staticdt=False, quiet=False):
        """Predict how a single colony will evolve in a single environment.

        Simulate the colony growing at a self-determined timestep up until
        one of the passed limters is reached.
        @param out string containing the name of the output file.
        @param tmax the maximum time for the simulation to run. Currently just over 1bn yr
        @param volmax the maximum volume fraction of the vessel the colony
          can take up.
        @param kmin the slowest rate of respiration to consider.
        @param endMF whether to close off the simulation when the average
          maintenance fraction in the colony exceeds 1.
        """

        resultsdict = self.dbh.buildresultsdict()
        logger.info('Checking if this is a new simulation')
        is_new_sim, self.SimID = self.dbh.createtable()

        if not is_new_sim:
            logger.info('This combination of organism(s) and location has been performed before. \n \n Aborting simulation, for results find: \n \t SimID: '+self.SimID)
            return


        logger.info('Setting up growth prediction.')
        stepmax = 100
        t=0
        if dt==None:
            dt = self.c.getmin_timestep(factorup=factorup)
            logger.debug('\n Simulation timestep selected: '+ str(dt))
        self.stoppingdict['Growth_Rate']['Min'] = self.stoppingdict['Growth_Rate']['Min']/dt


        # growth_watcher = [0]*len(self.c.collection)
        # growth_relpop = [0]*len(self.c.collection)
        # starting_gradient=[0]*len(self.c.collection)

        first=True
        timestart=timer.time()
        step=1
        while t < tmax:
            timestartstep=timer.time()
            step = step + 1#1+round(t/dt)

            logger.debug('\n Starting step '+ str(step))

            startpop = [o.get_population() for o in self.c.all()]# self.c.collection[0].num#self.c.get_alive_population()#*self.c.collection[0].volume
            try:
                self.c.take_step(dt)#, First=first) #
            except ZeroDivisionError as zde:
                if first:
                    logger.warning('Zero Division Error: Something may be out of food.')
                    logger.info('A Zero Division Error was encountered. Usually' + \
                      'this means that there is no power supply because the '+ \
                      'organism is trying to consume more than is available, ' + \
                      'maintenance fraction is above 1 in the first step, or ' + \
                      "reactions and reagents aren't defined properly")
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
                logger.info('Other error occurred while take step')
                logger.info('Saving data up to error...')
                for col in self.c.all():
                    col.output.appendvals(col.get_population(), dt)
                self.c.output.appendvals()
                resultsdict['Time'].append(t)
                resultsdict['Composition'].append(self.r.getconcs())
                full_results = self.dbh.dict_to_db(resultsdict, end=True)
                # print(full_results)
                raise e

            except e:
                raise e

            self.r.update_composition(dt)



            # Before, First did the current respiration pathway. Is it still needed?

            first=False
            # timetakestep=timer.time()
            # print("Taking step time: "+ str(timetakestep-timestartstep))

            t+=dt
            logger.info('Step '+ str(step) + ' complete\n')

            endpop = [o.get_population() for o in self.c.all()]


            if sum(endpop) > 0:
                # to format for sig figs use "{0:.6g}".format(data), 6 is no digits.
                # timestartdict=timer.time()
                resultsdict['Time'].append(t)
                resultsdict['Composition'].append(self.r.getconcs())

                # ake unique organism things the new dataheavy
                # if dataheavy:
                #     resultsdict = self.CRuniquemonitors(resultsdict, dt, startpop)


            else:
                full_results = self.dbh.dict_to_db(resultsdict, end=True)
                self.output_step_data(step, full_results)
                logger.info('All of the organisms have died')
                break


            if self.check_stopper('Maintenance_Fraction', self.c.get_maintenance_fractions(), resultsdict, step):
                break
            if self.check_stopper('Metabolic_Rate', self.c.get_metabolic_rates(), resultsdict, step):
                break
            if self.check_stopper('Growth_Rate', [abs(g) for g in self.c.get_growth_rates()], resultsdict, step):
                break
            if self.check_stopper('Volume_Fraction', [self.c.get_total_volume(inactive=True)], resultsdict, step):
                break
            if self.check_stopper('Population', [self.c.get_population(inactive=False)], resultsdict, step):
                break


            # if endMF and min(self.get_maintenance_fractions())>1.0:
            #     full_results = self.dbh.dict_to_db(resultsdict, end=True)
            #     self.output_step_data(step, full_results)
            #     logger.info('\n The average maintenance fraction has exceeded 1.')
            #     break
            #
            # ############### YOU ARE HERE ################
            #
            #
            # # will need updating for 2
            # if min(self.c.get_matabolic_rates())<kmin:
            #     full_results = self.dbh.dict_to_db(resultsdict, end=True)
            #     self.output_step_data(step, full_results)
            #     logger.info('The rate of respiration has slowed to below ' + \
            #     'your threshold for all organisms.')
            #     break
            #
            #
            #
            # if min(self.c.get_growth_rates())<=min_growthrate and step>200:
            #     # Either nothing is growing, or everything is dead, or growth
            #     # rate has dropped below your custom max growth rate
            #     full_results = self.dbh.dict_to_db(resultsdict, end=True)
            #     self.output_step_data(step, full_results)
            #     logger.info('The rate of growth has slowed to below ' + \
            #     'your threshold for all organisms.')
            #     break

            # but if none of these occur, output step data+save every 100 steps
            if (step/stepmax).is_integer():
                #update on where we are, and save data
                # print(resultsdict)

                logger.info('Saving Data')
                if t>=tmax:
                    full_results = self.dbh.dict_to_db(resultsdict, end=True)
                #elif step >=10000:
                    # this is a long simulation, so only save every 100th step
                #    full_results = self.dbh.dict_to_db(resultsdict, end=False, finonly=True)
                else:
                    full_results = self.dbh.dict_to_db(resultsdict, end=False)
                if not quiet:
                    self.output_step_data(step, full_results)

                resultsdict = self.dbh.buildresultsdict()

                if (step/10000).is_integer() and not staticdt:
                    # significant local changes are likely, adapt the timestep to compensate
                    dt = self.c.getmin_timestep(factorup=factorup)
                    logger.debug('\n Simulation timestep selected: '+ str(dt))
                    self.stoppingdict['Growth_Rate']['Min'] = self.stoppingdict['Growth_Rate']['Min']/dt



        logger.info('Simulation Complete')
        ex_time = timer.time()-timestart
        logger.info('This took '+str(ex_time)+' s  in real time')
        if not (step/stepmax).is_integer() and step != 1:
            logger.info('Saving Data')
            full_results = self.dbh.dict_to_db(resultsdict, end=True)
            if not quiet:
                self.output_step_data(step, full_results)


    def check_stopper(self, stopper_name, orglst, resultsdict, step):
        if min(orglst) > self.stoppingdict[stopper_name]['Max'] or max(orglst) < self.stoppingdict[stopper_name]['Min']:
            self.stoppingdict[stopper_name]['Count'] += 1
            if self.stoppingdict[stopper_name]['Count'] >= self.stoppingdict[stopper_name]['Consistency']:
                # full_results = self.dbh.dict_to_db(resultsdict, end=True)
                # self.output_step_data(step, full_results)
                logger.info('\n Stopped by ' + stopper_name+'.')
                return True
        else:
            self.stoppingdict[stopper_name]['Count'] = 0
            return False


    def output_step_data(self, step, dicti):#t, startpop, endpop):
        print('\n'+'step = ' + str(step))

        for key in sorted(dicti):
            print('---  '+key+'   '+str(dicti[key][-1]))


    #
    # def CRuniquemonitors(self, dicti, dt, startpop):
    #     TGP = [None]*len(self.c.collection)
    #     CO2W = [None]*len(self.c.collection)
    #     for i, OrgName in enumerate(self.dbh.OrgNames):
    #
    #         if OrgName == 'Methanogen':
    #             # add on the methanogen-specific monitors
    #             TGP[i] = self.c.collection[i].throttled_E_growth/(dt*startpop[i])
    #             # CO2W[i] = self.c.collection[i].respiration.rate*dt*((self.c.collection[i].P_growth-(self.c.collection[i].throttled_E_growth/dt))/(self.c.collection[i].P_s))
    #     # now put them into the dict as tuples
    #     for OrgName in self.dbh.OrgNames:
    #
    #         if OrgName == 'Methanogen':
    #             dicti['ThrottledGrowthPower'] = np.append(dicti['ThrottledGrowthPower'], tuple(TGP))
    #             dicti['WastedCO2'] = np.append(dicti['WastedCO2'], tuple(CO2W))
    #
    #     return dicti


    def empty(self):
        self.c.collection=[]
        self.r =self.backup_r
