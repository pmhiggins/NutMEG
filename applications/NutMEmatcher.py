from copy import deepcopy
import NutMEG as es

Tdef = ['Lever10pc', 'Lever2pc', 'Tijhuis']

class NutMEmatcher:
    """
    Class for matching parameters to outputs, specifically nutrients and maintenance.
    Cycle through implementations of ecosystem's predict_growth, until one gets
    the desired output.
    """

    def __init__(self, org, loc):
        self.org = org
        self.loc = loc

    def match(self, level='Basic',
      paramsdict={'MaintenanceFrac':[0,1]},
      target = {'GrowthRate':0.0003}):
        """Find an organism that matches the target output by varying paramsdict values.

        Parameters
        ----------
        level : str
            level of complexity to output. Basic: only efficiency. NutME: efficiency
            and the theoretical calculations which sould make it up.

        paramsdict : dict
            Which parameters to iterate through when matching. If a list with two
            entries is passed, iterate between those values. If a single value,
            stick to it, if more values use ONLY these values and return the
            best fit.
            TODO: options 2 and 3 here.

        target : dict
            Target paramter(s) to get to. At the moment this only supports
            growth rate, but others will be added.
            TODO: add other targets (final biomass? final conc?)

        """
        pd=deepcopy(paramsdict)
        t=deepcopy(target)

        e, SimID, OrgID = self.iterate_to_target(pd, t) #TODO: matching method which generates efficiency
        # send paramsdict, if NutrientFrac==[1], only do main, if not.....?

        if level == 'Basic':
            return e, SimID, OrgID
        elif level == 'NutME':
            self.print_theory_estimates(SimID, OrgID)
            return e, SimID, OrgID
        else:
            raise ValueError('Level unrecognised')


    def iterate_to_target(self, iterables, target):
        """Perform the iteration for the passed values.
        """
        for key, value in iterables.items():
            if len(value)==1:
                if key == 'MaintenanceFrac':
                    dt, cop = self.org.select_timestep(returncop=True)
                    initialPS = cop.P_s
                    self.org.maintenance.net_dict['Basal'] = value[0]*initialPS
                else:
                    raise ValueError('NutMEGmatch is not optimised for that parameter yet!')
            elif len(value)==2:
                if key == 'MaintenanceFrac':
                    return self.iterateMF(value, target)
                else:
                    raise ValueError('NutMEGmatch is not optimised for that parameter yet!')
            else:
                if key == 'MaintenanceFrac':
                    for i in value:
                        dt, cop = self.org.select_timestep(returncop=True)
                        initialPS = cop.P_s
                        self.org.maintenance.net_dict['Basal'] = value*initialPS
                else:
                    raise ValueError('NutMEGmatch is not optimised for that parameter yet!')


    def iterateMF(self, MFrange, target):
        """Iteration instructions for changing the maintenance fraction"""
        dt, cop = self.org.select_timestep(returncop=True)
        initialPS = cop.P_s
        met=False
        LID, OID, SID = None, None, None
        MF=1
        for k,v in target.items():
            while not met and MF>0.001:
                MF = (MFrange[1]+MFrange[0])/2.0

                org2 = deepcopy(self.org)
                loc2 = deepcopy(self.loc)
                org2.maintenance.net_dict['Basal'] = MF*initialPS

                LID, OID, SID = self.simulate_growth(org2, loc2)
                result = self.exparam(k, SID, OID)
                print(MF, result)
                hl = self.higherorlower(k, v, result)
                if hl =='+':
                    # need to increase MF for the next step.
                    MFrange[0] = MF
                elif hl == '-':
                    MFrange[1] = MF
                elif hl == ':-)':
                    MFrange = [MF,MF]
                    met=True
                print(MFrange)
        return MFrange[0], SID, OID



    def higherorlower(self, param, value, result):
        """For passed param, see if we have an acceptable result.
        Currently only works for GrowthRate and FinBM.
        """
        if 0.9 < result/value < 1.1:
            #regardless of value, if we're within 10%, call off optimisation
            return ':-)'
        if 'GrowthRate' or 'Volume' in param:
            # use the in syntax, in case something like GrowthRateVol is passed.
            # the growth rate should be approx independent of vol, cell, kg anyway.
            if result > value:
                #MF is TOO LOW.
                return '+'
            else:
                return '-'
        else:
            raise ValueError('NutMEGmatch is not optimised for that parameter yet!')


    def simulate_growth(self, org2, loc2):
        """Perfrom a growth prediction with the horde and locale passed"""
        org2.workoutID()
        loc2.dbh.workoutID()
        Cu = es.culture(hordes=[org2])

        ES = es.ecosystem(loc2, Cu)

        ES.predict_growth(quiet=True)

        return loc2.LocID, org2.OrgID, ES.dbh.SimID




    def print_theory_estimates(self, SimID, OrgID):
        """Print estimates at maintenance efficiency available in NutMEG."""
        PS, MF = self.net_powers(SimID, OrgID)
        org2 = deepcopy(self.org)
        print('Simulated Maintenance Fraction: '+str(MF))
        for Td in Tdef:
            org2.maintenance.Tdef = Td
            org2.maintenance.get_P_T()
            Tdv = org2.maintenance.net_dict['T']
            print('Theoretical T fraction: '+Td+': '+str(Tdv/PS))
        org2.maintenance.pHdef='FluxPerm'
        org2.maintenance.get_P_pH()
        print('Theoretical pH fraction: '+str(org2.maintenance.net_dict['pH']/PS))



    def net_powers(self, SimID, OrgID):
        """Return the power supply and maintenance fraction in the exponential
        phase of the simulation SimID."""
        PS = es.ecosystem_dbhelper.db_helper.extract_param_db_Sim(SimID, 'PowerSupply_'+OrgID)
        MF = es.ecosystem_dbhelper.db_helper.extract_param_db_Sim(SimID, 'MaintenanceFrac_'+OrgID)
        if 100 < 0.2*len(PS):
            return PS[100][0], MF[100][0]
        else:
            return PS[int(0.2*len(PS))][0], MF[int(0.2*len(PS))][0]

    def exparam(self, param, SimID, OrgID):
        """Return the requested simulation paramter in the exopnential phase of
        S=simulation SimID. """
        p = es.ecosystem_dbhelper.db_helper.extract_param_db_Sim(SimID, param+'_'+OrgID)
        if 'GrowthRate' in param:
            # to be taken in the exponential phase
            if 100 < 0.2*len(p):
                return p[100][0]
            else:
                return p[int(0.2*len(p))][0]
        elif 'Volume' in param:
            # to be taken at the end of the simulation
            return p[-1][0]
        else:
            raise ValueError('exparam does not recognise this parameter')
