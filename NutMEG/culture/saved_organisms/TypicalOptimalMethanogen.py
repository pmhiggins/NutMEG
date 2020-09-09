import NutMEG
import math, os
from statistics import mean
import NutMEG.util.NutMEGparams as nmp
from NutMEG.reactor.saved_systems.VenusDrop import VenusDrop
import pandas as pd
import numpy as np


#TODO: Tidy up
#TODO: add in base_organism


class TypicalOptimalMethanogen(NutMEG.base_organism):
    """Class for the typical optimal methanogen"""

    def __init__(self, R, orgtype='horde', paramchange={}, dbpath=nmp.std_dbpath, workoutID=False, **kwargs):
        params = TypicalOptimalMethanogen.avg_org_params(R, paramchange)
        params.update(kwargs)

        if orgtype == 'horde':
            NutMEG.horde.__init__(self,
              'TypicalOptimalMethanogen', R,
              self.setup_methanogenesis(R, params['k_RTP']),
              500,
              workoutID=workoutID,
              dbpath=dbpath,
              **params)


    def setup_methanogenesis(self, R, k_RTP=1.0):
        """return a methanogenesis reaction object in reactor R"""

        CO2 = NutMEG.reaction.reagent('CO2(aq)', R.env, phase='aq')
        H2aq = NutMEG.reaction.reagent('H2(aq)', R.env, phase='aq')
        CH4aq = NutMEG.reaction.reagent('CH4(g)', R.env, phase='g')
        H2O = NutMEG.reaction.reagent('H2O(l)', R.env, phase='l')

        thermalMG = NutMEG.reaction.reaction({CO2:1, H2aq:4}, {CH4aq:1, H2O:2},
          R.env)

        thermalMG.rate_constant_RTP = k_RTP

        return thermalMG


    @staticmethod
    def CH4Rate_from_GrowthRate(growthrate, CH4conc=1e-8):
        """Estimate the rate of CH4 production using the growth
        rate. ref Powell 1983."""
        return CH4conc*growthrate*(math.exp(growthrate) - 1)

    @staticmethod
    def avg_org_params(locale, paramchange={}, fromdata=False):
        """get the 'average' methanogen parameters from the methanogens csv."""
        # df = pd.read_csv(os.path.dirname(__file__)+'/../../data/methanogens.csv', header=0)
        # widths, lengths, Ts, pHs, CH4s, Pressures, GRs = [],[],[],[],[],[],[]
        # for index, row in df.iterrows():
        #     try:
        #         widths.append((float(row['Min. cell width'])+float(row['Max. cell width']))/2)
        #         lengths.append((float(row['Min. cell length'])+float(row['Max. cell length']))/2)
        #         # Ts.append(273.15+((float(row['Min. optimal growth temp.'])+float(row['Max. optimal growth temp.']))/2))
        #         # pHs.append((float(row['Min. optimal growth pH'])+float(row['Max. optimal growth pH']))/2)
        #         # Pressures.append(float(row['Pressure'])*1000)
        #         # GRs.append(row['Growth rate']/3600)
        #         # CH4s.append(TypicalOptimalMethanogen.CH4Rate_from_GrowthRate(row['Growth rate']/3600, CH4conc=paramchange.get('mol_CH4', 3e-8)))
        #     except Exception as e:
        #         print(str(e))
        #         continue

        # CH4s = np.polyfit(Ts, np.log(CH4s), 1)
        # GRs = np.polyfit(Ts, np.log(GRs), 1)
        #
        # GRT = math.exp(GRs[0]*(locale.env.T)+GRs[1])
        # CH4T = math.exp(CH4s[0]*(locale.env.T)+CH4s[1])

        # k_T = CH4T/(VenusDrop.getgasconc('CO2(aq)', 0.2*locale.env.P, locale.env.T, P_bar=locale.env.P, S=0)*(VenusDrop.getgasconc('H2(aq)', 0.8*locale.env.P, locale.env.T, P_bar=locale.env.P, S=0)**4))

        k_RTP = 0.
        maxmet = 0.
        if fromdata:
            TOMdf = pd.read_csv(os.path.dirname(__file__)+'/../../data/params2.csv')
            for index, row in TOMdf.iterrows():
                if TOMdf['Temperature'][index] == round(locale.env.T):
                    k_RTP = TOMdf['k_RTP'][index]
                    maxmet = TOMdf['MetabolicRate'][index]
        else:
            T = locale.env.T
            kpoly = [-5.31266542e-04,  4.46988203e-01, -9.67859408e+01]
            mpoly = [-4.07231645e-17,  1.32059394e-01, -8.10101975e+01]
            k_RTP = math.exp(kpoly[0]*(T*T)+kpoly[1]*(T)+kpoly[2])
            maxmet = math.exp(mpoly[0]*(T*T)+mpoly[1]*(T)+mpoly[2])

        vol=3.4412868852915668e-18
        # vol = (math.pi*((mean(widths)*0.5e-6)**2)*mean(lengths)*1e-6)
        # drym=vol*300
        # mass=vol*1000

        avg_org_uniqueparams = {'Tdef':paramchange.get('Tdef', 'Tijhuis'),
          'Basal':paramchange.get('Basal', 0),
          'volume':paramchange.get('volume', vol),
          'dry_mass':paramchange.get('dry_mass', vol*300),
          'mass':paramchange.get('mass', vol*1000), #'CH4rate':CH4T,
          'n_ATP':paramchange.get('n_ATP', 1.0),
          'lifespan':paramchange.get('lifespan',float('inf')),
          'k_RTP' : k_RTP,
          'max_metabolic_rate' : maxmet}

        return avg_org_uniqueparams
