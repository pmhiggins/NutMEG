

class Ecosystem:
    """
    Class containing an evolving population of residents (e.g., organisms),
    and their host chemical reactor. Behaves as the hub for time-integrated
    calculations or computing the state of a system containing multiple
    components.

    Attributes
    ----------
    reactor : Reactor
        Local chemical environment
    residents : list[Resident]
        List of evolving residents (e.g., different OrganismPopulation instances)
    """

    def __init__(self, reactor, residents):
        self.reactor = reactor
        self.residents = residents  # dict[str, Agent]

    def take_step(self, dt):
        """ Advance the ecosystem, its reactor and inhabitants by time dt in s. """
        for res in self.residents:
            res.take_step(dt, self)
        self.reactor.take_step(dt)
