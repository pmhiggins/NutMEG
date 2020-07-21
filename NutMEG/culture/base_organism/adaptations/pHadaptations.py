"""

Class for calculating adaptations to adverse pH conditions

"""
import math
# helper functions for Flux/permeability calcuation
def getfluxH(memb_pot, perm, pHint, pHext, area, K=38.7, T=298):
    return(perm*(T/298)*((K*memb_pot)/(1-math.exp(K*memb_pot)))*
      ((10**-pHext)-((10**-pHint)*math.exp(K*memb_pot)))*area*1e3)

def getfluxOH(memb_pot, perm, pHint, pHext, area, K=38.7, T=298):
    return(perm*(T/298)*((K*memb_pot)/(1-math.exp(K*memb_pot)))*
      ((10**(pHint-14))-((10**(pHext-14))*math.exp(K*memb_pot)))*area*1e3)

def getEnergyPump(pHint, pHext, T=298):
    return(8.31*T*abs(pHint-pHext))

class pHadaptations:


    def __init__(self, host):
        self.host = host # the host oragnism and its assosiated environment.

    """ Flux and Permeability calculation [QE Document] """


    def getFluxPerm_MP(self):
        # print(self.host.locale.composition['H+'].conc)
        pHext = -math.log10(self.host.locale.composition['H+'].conc)
        K= 96485/(8.31*self.host.locale.env.T)
        fluxH = getfluxH(self.host.memb_pot, self.host.PermH, self.host.pH_interior, pHext, self.host.surfacearea, K=K, T=self.host.locale.env.T)
        fluxOH = getfluxOH(self.host.memb_pot, self.host.PermOH, self.host.pH_interior, pHext, self.host.surfacearea, K=K, T=self.host.locale.env.T)
        E = getEnergyPump(self.host.pH_interior, pHext, T=self.host.locale.env.T)
        return abs(fluxH+fluxOH)*E
