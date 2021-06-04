"""

Class for calculating adaptations to adverse pH conditions

"""
import math
# helper functions for Flux/permeability calcuation
## out - in(exp)

## removed 1/T?

def GoldmanEQ(P, K, concs, SA, T, q=1):
    return P*(q)*(T/298)*(K/(1-math.exp(-q*K)))*(concs['in']-(concs['out']*math.exp(-q*K)))*SA*1e3
    # return P*(q**2)*(T/298)*(K/(1-math.exp(q*K)))*(concs['out']-(concs['in']*math.exp(q*K)))*SA*1e3

class pHadaptations:


    def __init__(self, host):
        self.host = host # the host oragnism and its assosiated environment.


    def _getfluxH(self):
        K = self.host.memb_pot*96485/(8.31*self.host.locale.env.T)
        concs = {'in':10**-self.host.pH_interior, 'out':10**-self.host.locale.pH}
        return GoldmanEQ(self.host.PermH, K, concs, self.host.surfacearea, self.host.locale.env.T, q=1)

    def _getfluxOH(self):
        K = self.host.memb_pot*96485/(8.31*self.host.locale.env.T)
        concs = {'in':10**(self.host.pH_interior-14),
          'out':10**(self.host.locale.pH-14)}
        return GoldmanEQ(self.host.PermOH, K, concs, self.host.surfacearea, self.host.locale.env.T, q=-1)

    def _getEnergyPump(self):
        return(8.31*self.host.locale.env.T*abs(self.host.pH_interior-self.host.locale.pH))


    """ Flux and Permeability calculation [QE Document] """

    def getFluxPerm_MP(self):

        fluxH = self._getfluxH()
        fluxOH = self._getfluxOH()
        E = self._getEnergyPump()
        return (abs(fluxOH)+abs(fluxH))*E
