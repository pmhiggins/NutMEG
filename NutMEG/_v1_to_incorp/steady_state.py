import NutMEG as es

from copy import deepcopy
from statistics import mean
import math, ast

class steady_state_1org:
    """
    Class for achieving a steady state within some bounding parameters.

    This application is a work-in-progress, so  not yet properly documented!
    """

    def __init__(self, org, loc, bounds,
      boundingparam='Quotient', timeframe=1e13):
        self.org = org
        self.loc = loc
        self.dbpath=self.org.dbh.dbpath

        # self.dim1 = {'vals': dim1, 'type':dim1type}
        # self.dim2 = {'vals': dim2, 'type':dim2type}
        self.timeframe = timeframe

        defaults = self.get_default_boundingparams()

        if boundingparam in defaults:
            self.bounds = bounds
            self.boundingparam = boundingparam
            self.boundingdict = defaults[boundingparam]
        else:
            raise ValueError('Unsupported boundingparam passed!')

    def get_default_boundingparams(self):

        return {'Quotient':{
          'orgloc':'loc', 'func':self.WithinQuotient,'stoppingdict':None},
          'Composition':{
            'orgloc':'loc', 'func':self.WithinComposition,'stoppingdict':{'Population':{'Max':float('inf'), 'Min':100, 'Consistency':100, 'Count':0}}}}

    def WithinQuotient(self, IDs):
        #Find out the quotient at the final step, and compare to the bounds
        if IDs['SimID'] == 'TooLong':
            return False
        try:
            DeltaG = es.ecosystem_dbhelper.db_helper.extract_param_db_Sim(
              IDs['SimID'],
              'EnergyAvailable_'+IDs['OrgID'],
              dbpath=self.org.dbh.dbpath)

            lnQ = (DeltaG[-1][0] - self.org.respiration.net_pathway.std_molar_gibbs)/(8.31*self.loc.env.T)
            if lnQ >= min(self.bounds) and lnQ <= max(self.bounds):
                return True
            else:
                return False
        except TypeError as e:
            print(e)
            print('')
            return False
        except Exception as e:
            print(e)
            return False

    def WithinComposition(self, IDs):
        # for this, bounds should be passed as a dict for the constaints you require.
        res = ''
        if IDs['SimID'] == 'TooLong':
            return 'TooLong'
        try:
            _Comp = es.ecosystem_dbhelper.db_helper.extract_param_db_Sim(IDs['SimID'], 'Composition', dbpath=self.org.dbh.dbpath)
            fincomp = ast.literal_eval(_Comp[-1][0])
            for key, value in self.bounds.items():
                if fincomp[key] >= min(value) and fincomp[key] <= max(value):
                    continue
                elif fincomp[key] <= 0.:
                    res += '!!'+key+' '
                else:
                    res += '!'+key + ' '
            return res
        except Exception as e:
            print(e)
            return 'Error! '+str(e)

    def is_ss(self, golden=False, dt=None):
        # return whether ss can be reached, and the IDs
        org2 = deepcopy(self.org)
        # loc2 = deepcopy(self.loc)
        IDs = self.simulate_growth(org2, org2.locale, golden=golden, dt=dt)
        ss_bool = self.boundingdict['func'](IDs)
        if ss_bool is True or type(ss_bool) == type(''):
            endMF = 0.
            endBM = 1.
            try:
                endMFa = es.ecosystem_dbhelper.db_helper.extract_param_db_Sim(IDs['SimID'], 'MaintenanceFrac_'+IDs['OrgID'], dbpath=self.org.dbh.dbpath)[-100:]
                endBM= es.ecosystem_dbhelper.db_helper.extract_param_db_Sim(IDs['SimID'], 'no_alive_'+IDs['OrgID'], dbpath=self.org.dbh.dbpath)[-1][0]
                endMF = mean([a[0] for a in endMFa])
            except:
                endMF = 1e10
            print(endMF, self.org.locale.ocean_pH,  self.org.locale.env.T)
            if endMF > 1. or endBM <= 100:
                return False, IDs
        return ss_bool, IDs

    def simulate_growth(self, org2, loc2, golden=False, dt=None):
        """Perfrom a growth prediction with the horde and locale passed"""

        org2.workoutID()
        loc2.dbh.workoutID()
        Cu = es.culture(hordes=[org2])

        ES = es.ecosystem(loc2, Cu, dbpath=self.dbpath)

        # global updates to stoppingdict to top brief stops.
        ES.stoppingdict['Maintenance_Fraction']['Max'] = float('inf')
        ES.stoppingdict['Growth_Rate']['Min'] = 0.
        ES.stoppingdict['Growth_Rate']['Consistency'] = 1000
        ES.stoppingdict['Volume_Fraction']['Max'] = float('inf')

        if self.boundingdict['stoppingdict']!=None:
            # there is a unique stoppingdict to use here!
            ES.stoppingdict.update(self.boundingdict['stoppingdict'])

        if dt== None:
            dt=1000.
            if org2.deathrate == 0.:
                dt=1000.
            else:
                dt= 0.01/org2.deathrate#max(1000000, 1/org2.deathrate)
                # dt = 0.01/org2.deathrate# Cu.getmin_timestep(factorup=1.01) # get timestep manually

        tmax = min([50000*dt, 1e9])#self.timeframe])
        if golden:
            dt=0.1/org2.deathrate
            tmax=900/org2.deathrate

        # tmax = 1e9#min([50000*dt, 1e7])#self.timeframe])

        if dt > tmax:
            return {'LocID':loc2.LocID, 'OrgID':org2.OrgID, 'SimID':'TooLong'}

        else:

            ES.predict_growth(quiet=True, dt=dt, tmax=tmax, staticdt=True)

            return {'LocID':loc2.LocID, 'OrgID':org2.OrgID, 'SimID':ES.dbh.SimID}
