import NutMEG
import pandas as pd
import os, math, sys, ast


class KineticallyLimitedOrganism(NutMEG.horde):

    def __init__(self, name, R, rxn,
      donor=None, acceptor=None,
      F_funcs=None,
      F_attrs=None,
      rate_func='zeroth order',
      n_ATP=1.0,
      num=500, **horde_kwargs):

        NutMEG.horde.__init__(self,
          name, R,
          rxn,
          num,
          n_ATP = n_ATP,
          **horde_kwargs)

        if donor:
            self.respiration.set_forcing_parameter('F_D', (lambda resp, K_D: resp.host.locale.composition[donor].conc/(resp.host.locale.composition[donor].conc + K_D)), ['K_D'])
        if acceptor:
            self.respiration.set_forcing_parameter('F_A', (lambda resp, K_A: resp.host.locale.composition[acceptor].conc/(resp.host.locale.composition[acceptor].conc + K_A)), ['K_A'])
        if F_funcs:
            for key, value in F_funcs.items():
                self.respiration.set_forcing_parameter(key, value[0], value[1])
        if F_attrs:
            for key, value in F_attrs.items():
                self.respiration.F_attrs[key] = value

        self.bm_conc = num*1e-12

        self.respiration.net_pathway.update_molar_gibbs_from_quotient(
          updatestdGibbs=False)

        self.update_metabolic_rate()

        self.donor_name = donor
        self.acceptor_name = acceptor


    @classmethod
    def Builtin(cls, ID, R, num=500, **horde_kwargs):

        orgs_df = pd.read_csv(os.path.dirname(__file__)+'/KLO_db.txt', sep='\t', index_col=0)

        if ID in orgs_df.index.tolist():#['Pathway'].tolist():


            rgts = {R.composition[k]:v for k,v in zip(ast.literal_eval(orgs_df.loc[ID]['ReactantNames']), ast.literal_eval(orgs_df.loc[ID]['ReactantStoichiometry']))}
            prods = {R.composition[k]:v for k,v in zip(ast.literal_eval(orgs_df['ProductNames'][ID]), ast.literal_eval(orgs_df['ProductStoichiometry'][ID]))}
            rxn = NutMEG.reaction.reaction(rgts, prods, R.env)
            rxn.rate_constant_env = orgs_df.loc[ID]['k_max']

            F_attrs = {k:orgs_df[k][ID] for k in ['K_A', 'K_D', 'xi']}

            F_funcs = {}
            if orgs_df['F_A_ID'][ID] != 'default':
                FFA, FFB = cls.builtin_forcing_funcs(orgs_df['F_A_ID'][ID], R)
                F_funcs['F_A'] = (FFA, FFB)
            if orgs_df['F_D_ID'][ID] != 'default':
                FFA, FFB = cls.builtin_forcing_funcs(orgs_df['F_D_ID'][ID], R)
                F_funcs['F_D'] = (FFA, FFB)


            return cls(ID, R, rxn,
              donor=orgs_df['Donor'][ID],
              acceptor=orgs_df['Acceptor'][ID],
              F_funcs=F_funcs,
              F_attrs=F_attrs,
              rate_func='zeroth order',
              n_ATP=orgs_df['n_ATP'][ID],
              num=num, **horde_kwargs)

        else:
            raise ValueError(ID+' not found in KLO_db.')

    @staticmethod
    def builtin_forcing_funcs(ID, R):

        if ID == 'MineralGoethite':
            return (lambda resp, K_A: 1./((resp.host.bm_conc/resp.host.locale.composition['Goethite'].conc) + K_A)), ['K_A']
        else:
            raise ValueError('Unknown custom forcing function bassed to builtin_forcing_funcs')

"""
class Methanogen_H2_HCO3(KineticallyLimitedOrganism):

    def __init__(self, R, name='Methanogen_H2_HCO3',
      F_funcs=None,
      F_attrs=None,
      n_ATP=0.125,
      **horde_kwargs):

          _F_attrs  = {'K_A':0., 'K_D':4.7e-6, 'xi':0.25}
          if F_attrs:
              _F_attrs.update(F_attrs) # overwrite / add as needed.

          rxn = NutMEG.reaction.reaction(
            {R.composition['H2(aq)']:1., R.composition['HCO3-']:0.25, R.composition['H+']:0.25},
            {R.composition['Methane(aq)']:0.25, R.composition['H2O(aq)']:0.75},
            R.env)

          KineticallyLimitedOrganism.__init__(self,
            name, R, rxn,
            donor='H2(aq)', acceptor='HCO3-',
            F_funcs=F_funcs,
            F_attrs=_F_attrs,
            n_ATP = n_ATP)


class Methanogen_Acetate(KineticallyLimitedOrganism):

    def __init__(self, R, name='Methanogen_Acetate',
      F_funcs=None,
      F_attrs=None,
      n_ATP=0.5,
      **horde_kwargs):

          _F_attrs  = {'K_A':0., 'K_D':2.3e-5, 'xi':1.0}
          if F_attrs:
              _F_attrs.update(F_attrs) # overwrite / add as needed.

          rxn = NutMEG.reaction.reaction(
            {R.composition['Acetate-']:1., R.composition['H2O(aq)']:1},
            {R.composition['HCO3-']:1, R.composition['Methane(aq)']:1},
            R.env)

          KineticallyLimitedOrganism.__init__(self,
            name, R, rxn,
            acceptor='Acetate-',
            F_funcs=F_funcs,
            F_attrs=_F_attrs,
            n_ATP = n_ATP)



class SulfateReducer_H2_SO4(KineticallyLimitedOrganism):

    def __init__(self, R, name='SulfateReducer_H2_SO4',
      F_funcs=None,
      F_attrs=None,
      n_ATP=0.25,
      **horde_kwargs):

          _F_attrs  = {'K_A':3.9e-5, 'K_D':1.1e-6, 'xi':1.5}
          if F_attrs:
              _F_attrs.update(F_attrs) # overwrite / add as needed.

          rxn = NutMEG.reaction.reaction(
            {R.composition['H2(aq)']:1., R.composition['SO4-2']:0.25, R.composition['H+']:0.25},
            {R.composition['HS-']:0.25, R.composition['H2O(aq)']:1.0},
            R.env)

          KineticallyLimitedOrganism.__init__(self,
            name, R, rxn,
            donor='H2(aq)', acceptor='SO4-2',
            F_funcs=F_funcs,
            F_attrs=_F_attrs,
            n_ATP = n_ATP)
"""
