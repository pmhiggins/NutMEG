from copy import deepcopy
import ast
import NutMEG as es
import NutMEG.util.NutMEGparams as nmp

class theory_estimates:
    """Class to easily work out some theory estimates for maintenance powers in
    specific environments

    Attributes
    ----------
    org : base_organism like
        The type of organism to use
    loc : reactor like
        The type of reactor (environment) to use
    """

    def __init__(self, org, loc):
        self.org = deepcopy(org)
        self.loc = deepcopy(loc)
        self.org.locale=self.loc
        self.dbpath=self.org.dbh.dbpath

    @classmethod
    def fromSim(cls, SimID, dbpath=nmp.std_dbpath):
        """Initialise from a Simulation already performed with ID SimID."""
        OrgIDs, LocID = es.db_helper.findOrgIDsLocID(SimID, dbpath=dbpath)
        #use the first organism
        OID = ast.literal_eval(OrgIDs)[0]

        return cls.fromOrgLoc(OID, LocID, dbpath=dbpath)

    @classmethod
    def fromOrgLoc(cls, OrgID, LocID, dbpath=nmp.std_dbpath):
        """Initialise from organisms or reactors already saved with IDs OrgID
        and LocID respectively."""
        R = es.reactor.r_from_db(es.db_helper.guess_name_from_ID(LocID)+'8020', LocID, dbpath=dbpath)
        Oname =  es.db_helper.guess_name_from_ID(OrgID)
        O = es.base_organism.bo_from_db(Oname, R, OrgID, dbpath=dbpath)

        return cls(O, R)



    def temperature_defenses(self, T, per_cell=True):
        """Return some expected temperature defense costs at temperature T which
        are built in to NutMEG, in units W/cell.

        Returns them in a dictionary, so far we have Lever10pc: Lever et al. (2015)
        ith protein replacement at 10% racemization, Lever2pc: the same with
        replacement at 2% racemization. and Tijhuis et al (1993)'s trend with
        empirical data.
        """
        ret = {'Lever10pc':0., 'Lever2pc':0., 'Tijhuis':0., 'TijhuisAerobe':0, 'TijhuisAnaerobe':0, 'Lever1/250':0}
        self.loc.change_T(T)
        for Td in ret.keys():
            self.org.maintenance.Tdef = Td
            self.org.maintenance.get_P_T()
            Tdv = self.org.maintenance.net_dict['T']
            if not per_cell:
                # return per unit volume biomass
                ret[Td] = Tdv/self.org.base_volume
            else:
                # return per cell
                ret[Td] = (Tdv)

        return ret

    #TODO: add in pH?
