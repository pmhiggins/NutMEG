import math
from .grower import Grower
from .maintainer import Maintainer
from .metaboliser import Metaboliser
from NutMEG.core.reactor.reaction import Reaction
from NutMEG.models.organism.growth_models.bioenergetic_growth_model import BioenergeticGrowthModel



class BaseOrganism:
    """
    This is a parent class for a model organism, to be used as the base
    reference point for organismic behaviour. It has a large number of attributes,
    but many are optional. The default is usually to behave like E. Coli.

    Attributes
    ----------
    name : str
        name of model organism, used for ecosystem management and output.
    metabolism : Metaboliser
        manage all metabolic behaviour such as rate of substrate uptake.
    maintenance : Maintainer, optional
        manage all maintnenace behaviour (Default is None, if None is passed on
        initiation the organism will have a maintainer with no survial costs)
    growth : Grower, optional
        Manages growth rates. Default behaviour is to use a Bioenergetic
        growth model.
    mass : float, optional
        Mass of each individual cell in kg. Default 1e-15
    dry_mass : float
        Dry mass of each individual cell in kg (e.g. its mass after being
        dehydrated) this is typically around 30% the mass. Default 3e-16
    volume : float
        volume of each individual cell in m^3. Default 1e-18
    age : float, kwarg
        Age of the organism in seconds. Default 0.
    E_synth : float, kwarg
        The amount of energy required to synthesise one cell in the locale.
        Default is 8x10^(-10) J (Higgins 2022). In NutMEG v1, it is possible
        to estimate using the synthesis module. For v2+ this is TBC.
    surfacearea : float, kwarg
        surface area of the cell im m^2 Default is to caluclate from volume
        assuming a spherical organism.
    base_life_span : float, kwarg
        Maximum age for the organism before it becomes inactive in s. Default
        is infinite.
    protein_fraction : float, kwarg
        Fraction of the organisms dry mass which is made up of proteins.
        Useful for synthesis calculations. Default is 0.55 (E. Coli)
    isactive : bool
        Whether the organism is active, e.g. interacting with the locale,
        metabolising, growing, etc. Starts as True, Changed to False when
        neccessary.
    """


    def __init__(self, name, metabolism,
      growth=None,
      maintenance=None,
      mass=1.e-15,
      dry_mass=3e-16,
      volume=1.e-18,
      surfacearea=None,
      E_synth=8e-10,#was 8e-10
      age=0.,
      base_life_span=None):
        """
        Parameters
        ----------
        name : str
            name of model organism, used for ecosystem management and output.
        metabolism : metaboliser
            manage all metabolic behaviour such as rate of substrate uptake.
        maintenance : Maintainer, optional
            manage all maintnenace behaviour (Default is None, if None is passed on
            initiation the organism will have a maintainer with no survial costs)
        growth : grower, optional
            Manages growth rates. Default behaviour is to use a Bioenergetic
            growth model.
        mass : float, optional
            Mass of each individual cell in kg. Default 1e-15
        dry_mass : float
            Dry mass of each individual cell in kg (e.g. its mass after being
            dehydrated) this is typically around 30% the mass. Default 3e-16
        volume : float
            volume of each individual cell in m^3. Default 1e-18
        age : float, kwarg
            Age of the organism in seconds. Default 0.
        E_synth : float, kwarg
            The amount of energy required to synthesise one cell in the locale.
            Default is 8x10^(-10) J (Higgins 2022). In NutMEG v1, it is possible
            to estimate using the synthesis module. For v2+ this is TBC.
        surfacearea : float, kwarg
            surface area of the cell im m^2 Default is to caluclate from volume
            assuming a spherical organism.
        base_life_span : float, kwarg
            Maximum age for the organism before it becomes inactive in s. Default
            is infinite.
        protein_fraction : float, kwarg
            Fraction of the organisms dry mass which is made up of proteins.
            Useful for synthesis calculations. Default is 0.55 (E. Coli)
        isactive : bool
            Whether the organism is active, e.g. interacting with the locale,
            metabolising, growing, etc. Starts as True, Changed to False when
            neccessary.
        """

        if name:
            self.name=name
            self.mass = mass
            self.dry_mass=dry_mass
            self.volume = volume
            self.surfacearea = surfacearea
            self.E_synth = E_synth #was 8e-10
            self.isactive=True
            self.issplitting=False
            self.base_life_span = base_life_span

        if not self.surfacearea:
            self.surfacearea = BaseOrganism.SAget(self.volume)

        if isinstance(metabolism, Metaboliser):
            self.metabolism = metabolism
        elif isinstance(metabolism, Reaction):
            self.metabolism = Metaboliser(metabolism)
        else:
            raise ValueError("base_organism cannot be initialised without a net metabolic reaction")

        if not maintenance:
            self.maintenance = Maintainer([])
        elif type(maintenance) is Maintainer:
            self.maintenance = maintenance
        else:
            warnings.warn('base_organism: '+name+' has been initalised without a maintainer')

        if not growth:
            self.growth = Grower()
        elif type(growth) is Grower:
            self.growth = growth
        else:
            warnings.warn('base_organism: '+name+' has been initalised without a grower')

    @classmethod
    def SAget(cls, v):
        """ estimate surface area from volume """
        return (4*math.pi*((v*3)**2))**(1/3)

    def update(self, locale):
        """
        Update the following time- and location-dependent organismic parameters:
        - metabolic rate
        - maintenance requirement (if present)
        - growth rate.
        This function does not update the reactor composition or free energies.
        """
        self.metabolism.compute_rate(self, locale)
        self.maintenance.compute_maintenance(self, locale)
        self.growth.compute_rate(self, locale)



    def get_metabolic_rate(self):
        """Return current metabolic rate (does not update it)."""
        return self.metabolism.rate

    def get_metabolic_equation(self):
        """Return current metabolic overall reaction equation."""
        return self.metabolism.net_pathway.equation

    def get_DeltaG(self, update=False, locale=None):
        """
        Return the molar Gibbs free energy of metabolism in J/mol

        Parameters
        ----------
        update : bool, optional
            Pass as True to update the free energy based on local context.
            Requires a locale is also passed.
        locale : Reactor, optional
            Local Reactor. Only required if update is True.
        """
        if update:
            return self.metabolism.update_DeltaG(locale)
        else:
            return self.metabolism.net_pathway.molar_gibbs


    def get_growth_rate(self):
        """Return current growth rate (does not update it)."""
        return self.growth.growth_rate


    def check_habitability(self, locale=None, update=True):
        """
        Perform bioenergetic habitability assessment for this organism.

        Returns a tuple containing: (bool of habitability result, Power supply, Maintenance power)

        Parameters
        ----------
        locale : Reactor, optional
            Local reactor, required if organism properties will be updated
        update : bool, optional
            Pass as True to update organism properties. Default True
        """
        if type(self.growth.growth_model) != type(BioenergeticGrowthModel()):
            raise ValueError('check_habtiability is not able to check non-bioenergetic habitability')

        if update:
            self.update(locale)

        p = self.growth.growth_model.cs_powers
        return p['P_s'] > p['P_m'], p['P_s'], p['P_m']


    # def get_ESynth(self, AA=False, comp=None):
    #     """Use the synthesis module to get the synthesis energy for this
    #     organism. Pass AA as True to include the cost of Amino Acid
    #     synthesis.
    #
    #     By default use E Coli parameters as built into synthesis. A future
    #     update might extend this, meaning we'll have to add comp """
    #
    #     #using comp in this way calls a  database which is already populated
    #     # passing host will only change results byt this organisms' protein_fraction
    #
    #     synth_dict = synth.get_ESynth_density(
    #       self.locale.env.T, AA=AA, compute=comp) #cost of sythesis per dry gram of cells
    #     #update E_synth for this organsim from the calculated density.
    #     if AA:
    #         self.E_synth = synth_dict*self.dry_mass*1000
    #     else:
    #         E_dens = []
    #         [E_dens.append(v) for k,v in synth_dict]
    #         self.E_synth = sum(E_dens)*self.dry_mass*1000

    def __str__(self):
        return self.name
