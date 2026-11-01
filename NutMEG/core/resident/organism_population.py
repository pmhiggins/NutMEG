from .resident import Resident

class OrganismPopulation(Resident):
    """
    Ecosystem resident for an evolving population of organisms.

    Attributes
    ----------
    name : str
        Identifier (e.g., 'Methanogen')
    base : BaseOrganism
        the model organism that will be used to compute bulk behavior
        of the organism population (e.g., metabolic and growth rates)
    num : float
        Number of active organisms
    death_rate : float
        Removal rate of organisms (unit s^-1), regardless of organismic properties
        like maintenance. Default 0.
    mass : float
        Total "wet" mass of this organism population in kg. Estimated using
        the base attribute.
    dry_mass : float
        Total dry mass of this organism population in kg. Estimated using
        the base attribute.
    """

    def __init__(self, base, num=1e3, death_rate=0.):
        self.age = 0.
        self.base = base
        self.num = num
        self.update_mass() # sets dry mass and "wet" mass using self.base

        self.death_rate = death_rate



    def take_step(self, dt, ES):
        """
        Advance the population forward in time by dt seconds in the local
        ecosystem ES. Grows/shinks the population, and performs the metabolic reaction in
        the reactor.

        Notes
        -----
        This default stepper only interacts with the
        ES.locale, but custom implementations could interact with other ES.residents


        Parameters
        ----------
        dt : float
            Time to advance the community in seconds.
        ES : Ecosystem
            Environmental context.
        """

        self.age += dt
        self.base.update(ES.reactor) # updates metabolic rate, maintenance + growth rate

        # perform the catabolic reaction with the locale
        moles_consumed = self.num*self.base.get_metabolic_rate()*dt
        ES.reactor.perform_reaction(self.base.get_metabolic_equation(), moles_consumed)
        for k,v in self.base.metabolism.extra_pathways:
            ES.reactor.perform_reaction(v[0], self.num*v[1]*dt)

        # update num based on growth
        net_gr = max(-1., self.base.get_growth_rate() - self.death_rate)
        new_cells = self.num * net_gr * dt
        self.num += new_cells
        self.update_mass()


    def update_mass(self):
        """ Update mass and dry_mass using cell-specific data in ``base``."""
        self.mass = self.num * self.base.mass
        self.dry_mass = self.num * self.base.dry_mass

    def get_num(self):
        return self.num
