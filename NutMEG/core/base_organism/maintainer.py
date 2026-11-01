
class Maintainer:
    """
    This class is for computing and calculating the maintenance
    requirements in the form of powers for a given organism.

    All values are per cell.

    Attributes
    ----------
    mechanisms : list[MaintenanceModel]
        List of MaintenanceModel object representing the different maintenance
        stresses the cell faces.
    total_maintenance_power : float
        Aggregate sum of all maintennce power requirements for host cell.
    MP_list : list
        List same length as ``mechanisms`` detailing the indidual contributions
        to total maintenance of each mechanism modelled.
    rebuild : list, optional
        YET TO BE IMPLEMENTED IN v2.
        Which maintenance processes required biomass replacement. List should
        contain string keys corresponding to the maintenance processes that
        require biomass synthesis and hence reduce the energy available for
        new growth. For example, setting ``rebuild = ['T']`` tells the maintainer
        that the energy used to defend against temperature must be used to
        synthesise biomass and hence also consume the organism's nutrient budget.
        Default value is empty (ie, maintenance processes do not need new biomass).
    """

    #P_store = 0. # Fractional power stored away for future use. Depreciated,
    # but I'll reintroduce it if bugs appear.


    def __init__(self,
      mechanisms=[]):
      # rebuild=[]):
        """
        Parameters
        ----------
        mechanisms : list[MaintenanceModel]
            List of MaintenanceModel object representing the different maintenance
            stresses the cell faces.
        """
        self.mechanisms = mechanisms

        self.total_maintenance_power = 0.
        self.total_rebuild_power = 0.

        self.MP_list = [0,]*len(mechanisms)

        # self.rebuild = rebuild



    def add_mechanism(self, MM):
        """Add a mechanism to the list of MaintenanceModels"""
        self.mechanisms.append(MM)
        self.MP_list.append(0.)

    def compute_maintenance(self, host, locale):
        _PM = 0.
        for i, MM in enumerate(self.mechanisms):
            this_PM = MM.compute(host,locale)
            self.MP_list[i] = this_PM
            _PM += this_PM
        self.total_maintenance_power = _PM
        return _PM


    # def get_P_store(self):
    #     """Get the net power stored if there is any"""
    #     if 'store' in self.net_dict:
    #         return self.net_dict['store']
    #     else:
    #         return 0.0
    #
    # def calculate_P_rebuild(self):
    #     """
    #     sum together all the maintenance powers that are flagged as
    #     requiring biomass rebuilt
    #     """
    #
    #     self.P_rebuild = 0.
    #     for rb in self.rebuild:
    #         self.P_rebuild += self.net_dict[rb]
