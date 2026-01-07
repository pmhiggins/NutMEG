import NutMEG
import pandas as pd
import os, math, sys, ast, yaml


class KineticallyLimitedOrganism(NutMEG.horde):
    """
    Class for establishing an organism with known kinetic forcing parameters.

    The KLO initialises a kinietically limited horde-like object in a more
    user-friendly way than the base horde class. Forcing functions and their
    attributes can be passed directly to the __init__ function or, more simply,
    use the Builtin constructor function to read-in a known
    KineticallyLimitedOrganism from one of NutMEG's kinetic databases.

    # TODO: add explainer for kinetic databases to NutMEG docs.

    """

    def __init__(self, name, R, rxn, num,
      kinetic_F_funcs=None,
      kinetic_F_attrs=None,
      growth_F_funcs=None,
      growth_F_attrs=None,
      donor=None, acceptor=None,
      n_ATP=1.,
      kinetic_rate_func='zeroth order',
      rate_constant_env=None,
      max_growth_rate=None,
      **horde_kwargs):

        print(horde_kwargs)

        respiration_kwargs={
          'kinetic_forcing_parameters': kinetic_F_funcs,
          'kinetic_F_attrs': kinetic_F_attrs,
          'rate_constant_env': rate_constant_env,
          'rate_func':kinetic_rate_func,
          'n_ATP':n_ATP}

        horde_kwargs['respiration_kwargs'] = respiration_kwargs

        CHNOPS_kwargs = {
          'max_growth_rate' : max_growth_rate,
          'CHNOPS_forcing_parameters': growth_F_funcs,
          'CHNOPS_F_attrs': growth_F_attrs,
        }

        horde_kwargs['CHNOPS_kwargs'] = CHNOPS_kwargs

        super().__init__(
          name, R,
          rxn,
          num,
          **horde_kwargs)

        # # if donor:
        # #     self.respiration.set_forcing_parameter('F_D', (lambda resp, K_D: resp.host.locale.composition[donor].conc/(resp.host.locale.composition[donor].conc + K_D)), ['K_D'])
        # # if acceptor:
        # #     self.respiration.set_forcing_parameter('F_A', (lambda resp, K_A: resp.host.locale.composition[acceptor].conc/(resp.host.locale.composition[acceptor].conc + K_A)), ['K_A'])
        #
        # # should just pass these to the respirator constructor!
        # if kinetic_F_funcs:
        #     for key, value in kinetic_F_funcs.items():
        #         self.respiration.set_forcing_parameter(key, value[0], value[1])
        # if kinetic_F_attrs:
        #     for key, value in kinetic_F_attrs.items():
        #         self.respiration.F_attrs[key] = value
        #
        # if max_growth_rate:
        #     self.max_growth_rate = max_growth_rate
        #     # unless explicitly specified otherwise, tell growth calcs to
        #     # consider nutrient limitation
        #     self.nutrientlimitation = horde_kwargs.get('nutrientlimitation', True)
        # if growth_F_funcs:
        #     for key, value in growth_F_funcs.items():
        #         self.CHNOPS.set_forcing_parameter(key, value[0], value[1])
        # if growth_F_attrs:
        #     for key, value in growth_F_attrs.items():
        #         self.CHNOPS.F_attrs[key] = value

        self.bm_conc = num*self.mass*1000 # g / L

        # # not needed if building in the respirator constructor
        # self.respiration.net_pathway.update_molar_gibbs_from_quotient(
        #   updatestdGibbs=False)
        # self.respiration.get_rate()

        self.donor_name = donor
        self.acceptor_name = acceptor


    @classmethod
    def Builtin(cls, ID, R, num=500, db_fn='default', **horde_kwargs):

        print(ID)
        org_date = None

        if db_fn == default:
            db_fn = os.path.dirname(__file__)+'/KLO_db.yaml'

        with open(db_fn, 'r') as f:
            org_data = yaml.safe_load(f)

        org_props = None
        try:
            org_props = org_data.get(ID)
        except:
            raise ValueError(ID+' not found in KLO_db.')


        rgts = {R.composition[k]:v for k,v in org_props.get('Reactants', {}).items()}
        prods = {R.composition[k]:v for k,v in org_props.get('Products', {}).items()}

        rxn = NutMEG.reaction.reaction(rgts, prods, R.env)

        rxn.rate_constant_env = org_props.get('Bioenergetics')['k_max'] * org_props.get('dry_mass', 3e-13)

        # forcing function attributes
        kinetic_F_attrs = {
          'xi' : org_props.get('Bioenergetics', {}).get('xi')
        }


        # look for any additional kinetic forcing functions
        kinetic_F_funcs = {}
        # org_kinetic_forcing = org_props['KineticForcing'] # a dict containing kinetic forcing data

        # potential error source below. If doesn't work go back to commented out above.
        # Trying to make this work when there is no KineticForcing of GrowthForcing
        # so the KLO can be more flexible
        org_kinetic_forcing = org_props.get('KineticForcing', {})

        for F_k, F_vs in org_kinetic_forcing.items():

            this_func_attrs = {}
            for k, v in F_vs.items():
                # this is a dict of information relevant to this forcing parameter,
                # identifying it, citing a source, and providing the
                # relevant attributes (e.g., for a Monod model: S and K)
                if k not in ['F_ID', 'func', 'ref', 'note', 'notes']:
                    # define a unique identifier for this attribute
                    _ID = k+'_'+F_vs['F_ID']+F_vs['func']

                    # assign its value to F_attrs, which will be stored in the organism
                    kinetic_F_attrs[_ID] = v

                    # connect the generic attribute names to the unique ones
                    # e.g., for monod: 'S':'S_DonorMonod' etc
                    this_func_attrs[k] = _ID

            # setup the forcing function for this mechanism
            FFA, FFB = cls.builtin_forcing_funcs(F_vs['func'], this_func_attrs)
            kinetic_F_funcs[F_k] = (FFA, FFB)

            # F_ID = org_props.get(F, {'ID':'default'}).get('ID')
            # # {'ID':'default'} is used as the catch to make sure changes
            # # aren't made where the forcing function is not applicable.
            #
            # if F_ID != 'default':
            #     FFA, FFB = cls.builtin_forcing_funcs(F_ID, R)
            #     F_funcs[F] = (FFA, FFB)

            # repeat the process above for growth forcing functions
            growth_F_funcs = {}
            growth_F_attrs = {}
            org_growth_forcing = org_props.get('GrowthForcing', {}) # a dict containing growth forcing data
            for F_k, F_vs in org_growth_forcing.items():
                this_func_attrs = {}
                for k, v in F_vs.items():
                    if k not in ['F_ID', 'func', 'ref', 'note', 'notes']:
                        _ID = k+'_'+F_vs['F_ID']+F_vs['func']
                        growth_F_attrs[_ID] = v
                        this_func_attrs[k] = _ID
                FFA, FFB = cls.builtin_forcing_funcs(F_vs['func'], this_func_attrs)
                growth_F_funcs[F_k] = (FFA, FFB)


        return cls(ID, R, rxn, num,
          donor=org_props.get('Donor'),
          acceptor=org_props.get('Acceptor'),
          kinetic_F_funcs=kinetic_F_funcs,
          kinetic_F_attrs=kinetic_F_attrs,
          growth_F_funcs=growth_F_funcs,
          growth_F_attrs=growth_F_attrs,
          kinetic_rate_func='zeroth order',
          n_ATP=org_props.get('Bioenergetics').get('n_ATP', 1.0),
          max_growth_rate=org_props.get('Bioenergetics').get('mu_max', None),
           **horde_kwargs)



    # @staticmethod
    # def builtin_forcing_funcs(funcID, attrs):
    #
    #     if funcID == 'Monod':
    #         # 2 attrs: substrate ID (e.g., 'H2(aq)'), and Monod half-saturation constant
    #         return (lambda resp, S, K: resp.host.locale.composition[S].conc/(resp.host.locale.composition[S].conc + K), [attrs['S'], attrs['K']])
    #     if funcID == 'MineralGoethite':
    #         return (lambda resp, K: (resp.host.bm_conc/resp.host.locale.composition['Goethite'].conc)/((resp.host.bm_conc/resp.host.locale.composition['Goethite'].conc) + K)), [attrs['K']]
    #
    #     else:
    #         raise ValueError('Unknown custom forcing function bassed to builtin_forcing_funcs')

    @staticmethod
    def suggest_orgs(R, products=False):
        """
        Suggest a list of organism keys in the KLO database that are
        compatible with reactor R for habitability analyses.
        """
        with open(os.path.dirname(__file__)+'/KLO_db.yaml', 'r') as f:
            data = yaml.safe_load(f)
        print(data)

        viables = []
        for org in data.keys():
            viable = True
            for k,v in data[org]['Reactants'].items():
                # k is the species name
                if not R.contains_reagent(k):
                    viable=False
            if products:
                for k,v in data[org]['Products'].items():
                    # k is the species name
                    if not R.contains_reagent(k):
                        viable=False
            if viable:
                viables.append(org)

        return viables


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
