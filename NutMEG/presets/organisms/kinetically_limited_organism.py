
from ... import core as nmc
from ...core import Reaction as rxn

from ...models.organism.forcing_factors import find_ff
from ...models.organism.base_rate_models import find_brm
from ...models.organism.growth_models import find_gm

import yaml, os


class KineticallyLimitedOrganism:
    """
    Class for creating BaseOrganism objects using databases of pre-defined
    organism specific kinetic and thermodynamic properties. 
    """

    builtin_KLO_dbs = ['example_KLO_db.yaml', 'Dale2006.yaml']

    def get(R, org_key, KLO_db='example_KLO_db.yaml', adjustments={}):
        """
        Return a BaseOrganism object constructued using a database of
        organism-specific kinetic and thermodynamic parameters.

        Parameters
        ----------
        R : Reactor
            Local reactor to initialise the organism (and its metabolism) in.
        org_key : str
            String identifier for the organism in the database
        KLO_db : str, optional
            Name of built-in database to use, or the path to a user-defined one.
            The database must be in yaml format.
        adjustments : dict, optional
            A dictionary which will be used to overwrite organism parameters as
            they are in the database and/or upon initialisation.
        """

        # retrieve the KLO database and assign dict to org_data
        KLO_fn = ''
        org_data = None

        if KLO_db in KineticallyLimitedOrganism.builtin_KLO_dbs:
            KLO_fn = os.path.dirname(__file__)+'/../../data/presets/KLO/'+KLO_db
        else:
            KLO_fn = KLO_db

        with open(KLO_fn, 'r') as f:
            org_data = yaml.safe_load(f)

        # retrieve entry for target organism
        org_props = None
        try:
            org_props = org_data.get(org_key)
        except:
            raise ValueError(org_key+' not found in '+KLO_db)

        # add adjustments to database organism as passed by user.
        org_props = KineticallyLimitedOrganism.deep_update(org_props, adjustments)

        # set up net metabolic reaction
        rgts = org_props.get('Reactants', {})
        prods = org_props.get('Products', {})

        met_rxn = rxn(R, rgts, prods, add_missing_rgts=True)

        # Read in parameters for the Metaboliser

        Met_dict = org_props.get('Metabolic', {})

        kinetic_brm = 'default'
        if 'Rate' in Met_dict:
            _brm = find_brm(Met_dict['Rate']['func'])
            kinetic_brm = _brm(**Met_dict['Rate'])

        kinetic_ff_objs = []
        kinetic_ff_labels = []
        for k,v in Met_dict.get('Forcing', {}).items():
            kinetic_ff_labels.append(k)
            _ff = find_ff(v['func'])
            kinetic_ff_objs.append(_ff(**v))

        _met = nmc.org.Metaboliser(met_rxn,
          base_rate = kinetic_brm,
          forcing_factors = kinetic_ff_objs,
          forcing_factor_labels = kinetic_ff_labels)

        # Read in parameters for the Grower
        Gro_dict = org_props.get('Growth', {})


        growth_brm = 'default'
        if 'Rate' in Gro_dict:
            _brm = find_brm(Gro_dict['Rate']['func'])
            growth_brm = _brm(**Gro_dict['Rate'])

        growth_gm = 'default'
        if 'Model' in Gro_dict:
            _gm = find_gm(Gro_dict['Model']['func'])
            growth_gm = _gm(**Gro_dict['Model'])

        # Read in forcing_factors for the Grower (if any)
        growth_ff_objs = []
        growth_ff_labels = []
        for k,v in org_props.get('Growth', {}).items():
            growth_ff_labels.append(k)
            _cls = find_ff(v['func'])
            growth_ff_objs.append(_cls(**v))

        _gro = nmc.org.Grower(base_rate=growth_brm,
          forcing_factors = growth_ff_objs,
          forcing_factor_labels = growth_ff_labels,
          growth_model = growth_gm)

        return nmc.BaseOrganism(org_key, _met, _gro, **org_props.get('BaseOrganism', {}))




    @staticmethod
    def deep_update(d, u):
        """
        Update dictionary d with dictionary u, mindful of depth. e.g., this will
        update nested dictionary keys without affecting other keys.
        """
        for ak, av in u.items():
            if isinstance(av, collections.abc.Mapping):
                d[ak] = KineticallyLimitedOrganism.deep_update(d.get(ak, {}), av)
            else:
                d[ak] = av
        return d
