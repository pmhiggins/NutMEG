from .GrowthModel import GrowthModel

class BioenergeticGrowthModel(GrowthModel):
    """
    Superclass for all growth rate calculations.

    Attributes
    ----------
    requires : list
        List of additional required host properties to run
        this GrowthModel.
    """

    def __init__(self, host, locale):
        """
        Parameters
        ----------
        host : ``base_organism'' like
            Host organism. Some GrowthModels will need this to initialise and
            some won't. It is best to assume they will (else they may throw an error)
        locale : ``reactor'' like
            Host chemical reactor. Some GrowthModels will need this to initialise and
            some won't. It is best to assume they will (else they may throw an error)
        """
        super().__init__(host,locale)
        self.requires = {'G_C':host.metabolism.forcing_factors}
        self.cs_powers = {'P_s': None, 'P_m':None, 'P_g_gross':None, 'P_g_net':None}


    def compute(self, host, locale):
        """ Calculate and return the growth rate according to this model."""

        G_C = self.find_subparams()['G_C']
        try:
            host.growth.max_growth_rate = host.growth.base_rate.compute(host, locale)
        except NotImplementedError:
            host.growth.max_growth_rate = float('inf')

        values = [f.compute(host, locale) for f in host.growth.forcing_factors]
        expected_growth_rate =  host.growth.aggregator.combine(host.growth.max_growth_rate, values)
        expected_P_growth = host.E_synth * expected_growth_rate

        P_growth = (G_C * host.metabolism.rate) - host.maintenance.total_maintenance_power
        P_rebuild = host.maintenance.total_rebuild_power

        net_growth_P = None
        gross_growth_P = None

        if (P_rebuild + P_growth) < expected_P_growth:
            # the biomass growth is energy limited
            net_growth_P = P_growth #expected_P_growth - P_rebuild
            gross_growth_P = P_rebuild + P_growth
        else:
            # the energy-limited biomass cannot be built in this step. Other
            # forcing factors are controlling it.
            # instead, build the rate-limited amount P_ex.
            if P_rebuild <= expected_P_growth:
                # maintenance is possible, but not the originally intended growth.
                net_growth_P = expected_P_growth - P_rebuild
                gross_growth_P = expected_P_growth
            else:
                # not even complete maintenance is possible, a negative net growth
                # power is needed (i.e. the biosphere will shrink)
                net_growth_P = expected_P_growth - P_rebuild
                gross_growth_P = expected_P_growth

            # reduce the metabolic rate to accord to the actual
            # power supply that is useable.
            host.metabolism.rate -= (P_growth - expected_P_growth)/G_C

        host.growth.growth_rate = net_growth_P / host.E_synth
        host.growth.gross_growth_rate = gross_growth_P / host.E_synth

        self.cs_powers = {
            'P_s' : host.metabolism.rate * G_C,
            'P_m' : host.maintenance.total_maintenance_power,
            'P_g_gross' : gross_growth_P,
            'P_g_net' : net_growth_P,
        }

        return host.growth.growth_rate, host.growth.gross_growth_rate

    def outputs(self):
        return {'cs_powers': self.cs_powers}
