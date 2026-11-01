from ..core import Reagent as rgt
from ..core import Reactor as rtr
from ..core import Reaction as rxn


import sys, os
import pandas as pd
from uncertainties import ufloat as uf


class Enceladus:
    """
    Class for building Reactors that resemble Enceladus's subsurface ocean.
    """

    Waite2017ratios_uf = {'CO2': uf(0.55, 0.25), 'CH4': uf(0.2, 0.1), 'NH3': uf(0.85, 0.45), 'H2': uf(0.9,0.5), 'H2S':uf(0.0021,0.001)}
    Waite2017ratios_nom = {'CO2': 0.55, 'CH4': 0.2, 'NH3': 0.85, 'H2': 0.9, 'H2S':0.0021}

    H24_species = ['H2O(aq)','H+','OH-','Na+','Cl-','HCO3-','CO2(aq)','CO3-2']

    def get_Enceladus(key, enc_kwargs={}, rtr_kwargs={}):
        """
        Return a Reactor object mimicking an Enceladus environment. Multiple
        models will be available, and the model of choice can be selected using the
        `key` parameter. Options currently include:

        'Higgins2024' : The chemical speciations reported in Higgins et al (2024)
        JGR: Planets. Available at temperatures between 273-473 K, pressures 1bar
        and 100bar, and [Cl] 0.05 to 0.4 (Corresponding to [DIC] also, see the
        paper for details).

        Parameters
        ----------
        key : str
            Key determining which model to build the Enceladus Reactor from.
        enc_kwargs : dict, optional
            kwargs to pass to the specific Enceladus model
        rtr_kwargs : dict, optional
            kwargs to pass to the Reactor initialisation.
        """

        if key == 'Higgins2024':
            return Enceladus.get_Enceladus_Higgins2024(enc_kwargs, **rtr_kwargs)


    def get_Enceladus_Higgins2024(pH_bo=8., Cl=0.1, P=1,
      spec_model='pitzerPHREEQC',
      mixingratios='default',
      chemical_species = 'default',
      from_fn = None,
      **rtr_kwargs):
        """
        Return a Reactor containing an Enceladus-like composition based on
        those calculated in Higgins et al., (2024) JGR:Planets.

        Only a subset of the chemical speciations in that work are available.
        To read in additional speciation results, they can be downloaded from
        the Higgins et al., (2024) dataset and the filename passed to the kwarg
        `from_fn`.

        Parameters
        ----------
        pH_bo : float, optional
            pH of the ocean at 273.15 K. Default 8.
        Cl : float, optional
            Molality of Cl, which defines the chemical speciation output to
            retrieve. Values can be either 0.05, 0.1, 0.2, 0.4. For more fidelity
            in Cl molality, manually read in chemical speciation data from the
            zenodo archive referenced in Higgins et al., (2024).
        P : float, optional
            Pressure in bar. Value can be 1 or 100 to match with the Higgins
            et al (2024) speciation endmembers.
        mixingratios : dict, optional
            Dict of gas mixing ratios to incorporate into the ocean. Default
            is to use the nominal values of the mixing ratios reported by
            Waite et al 2017 Science.
        chemical_species : list, optional
            List of chemical species to extract from the Higgins et al 2024
            chemical speciation and add to the returned ``Reactor``.
        from_fn : str, optional
            Path to filename of a speciation output from Higgins et al. (2024).
            If this is not None, that file will be used regardless of other Cl
            or P values passed. pH_bo (and T if not 273K) is still required.
        **rtr_kwargs : dict
            Additional kwargs for the ``Reactor`` being built. Most relevant
            to this preset is the temperature 'T'.
        """

        _mixingratios = None
        if mixingratios == 'default':
            _mixingratios = Enceladus.Waite2017ratios_nom
        elif type(mixingratios) is type({}):
            _mixingratios = mixingratios
        else:
            raise ValueError('Unrecognised form of mixingratios passed to get_Enceladus_Higgins2024')

        if chemical_species == 'default':
            chemical_species = Enceladus.H24_species

        EncDefaults = {
          'thermodb' : 'supcrt07-organics', # for max replicability with H+24. SUPCRT was used as no CH4 in phreeqc
          'T' : 273.15,
          'P' : P*1e5,
          'kgH2O' : 1.}
        EncDefaults.update(rtr_kwargs)

        _this = rtr(**EncDefaults)


        ### retrieve the chemical speciation
        spec_fn = None
        if from_fn:
            # use user-defined spec file
            spec_fn = from_fn
        else:
            # use builtin spec file
            specdir = os.path.dirname(__file__)+'/../data/presets/Enceladus/H24speciation/Clconc_'+str(Cl)
            spec_fn = specdir+'/spec_'+str(int(P))+'bar_'+spec_model+'.csv'
        df = pd.read_csv(spec_fn)
        _df = df[df['T'] == _this.T]
        __df = _df[_df['pH_bo'] == pH_bo] #isolates for the line of the df we need.


        for v in chemical_species:
            _v = v
            if v == 'H2O(aq)':
                _v = 'H2O'
            # the species is not present in the reactor yet
            g = float(__df['g'+_v].iloc[0])
            m = float(__df['m'+_v].iloc[0])
            rct = rgt(v, _this, phase='aq',
              amount=(m,'molal'), activity=m*g, gamma_molal=g)


        ### add the dissolved gases H2 and CH4
        mol_CH4 = (_mixingratios['CH4']/_mixingratios['CO2'])*_this.composition['CO2(aq)'].get_molality()

        CH4aq = rgt('Methane(aq)', _this, phase='aq', amount=(mol_CH4, 'molal'))

        mol_H2 = (_mixingratios['H2']/_mixingratios['CO2'])*_this.composition['CO2(aq)'].get_molality()

        H2aq = rgt('H2(aq)', _this, phase='aq', amount=(mol_H2, 'molal'))

        _this.update_pH(_from='H+')
        return _this


    def add_methanogenesis(enc_reactor):
        """
        Add methanogenesis in the form 4H2 + CO2 -> CH4 + 2H2O to reactor
        ``enc_reactor``.
        """

        r = {'CO2(aq)':1, 'H2(aq)':4}
        p = {'Methane(aq)':1, 'H2O(aq)':2}
        MG = rxn(enc_reactor, r, p)
        return MG
