from copy import deepcopy
import ast
import NutMEG as es
import NutMEG.util.NutMEGparams as nmp

class theory_estimates:

    def __init__(self, org, loc):
        self.org = deepcopy(org)
        self.loc = deepcopy(loc)
        self.org.locale=self.loc
        self.dbpath=self.org.dbh.dbpath

    @classmethod
    def fromSim(cls, SimID, dbpath=nmp.std_dbpath):
        OrgIDs, LocID = es.db_helper.findOrgIDsLocID(SimID, dbpath=dbpath)
        #use the first organism
        OID = ast.literal_eval(OrgIDs)[0]

        return cls.fromOrgLoc(OID, LocID, dbpath=dbpath)
        """
        R = es.reactor.r_from_db(es.db_helper.guess_name_from_ID(LocID), LocID, dbpath=dbpath)

        Oname =  es.db_helper.guess_name_from_ID(OID)
        O = es.base_organism.bo_from_db(Oname, R, OID, dbpath=dbpath)
        return cls(O, L)
        """

    @classmethod
    def fromOrgLoc(cls, OrgID, LocID, dbpath=nmp.std_dbpath):
        R = es.reactor.r_from_db(es.db_helper.guess_name_from_ID(LocID)+'8020', LocID, dbpath=dbpath)
        Oname =  es.db_helper.guess_name_from_ID(OrgID)
        O = es.base_organism.bo_from_db(Oname, R, OrgID, dbpath=dbpath)

        return cls(O, R)



    def temperature_defenses(self, T, per_cell=True):
        ret = {'Lever10pc':0., 'Lever2pc':0., 'Tijhuis':0.}

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
