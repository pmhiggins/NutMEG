"""

Class for calculating adaptations to adverse T conditions

"""
import math

class Tadaptations:


    def __init__(self, host):
        self.host = host # the host oragnism and its assosiated environment.

    """ Lever calculation [QE Document] """

    def getLeverME(self, cutoff_pc=10):
        """Power cost due to temperature according to Lever et al 2015
        @param cutoff_pc the % racemization of AA in proteins at which it is
          replaced."""

        k_yr = 0.00012*math.exp(0.10174*(self.host.locale.env.T-273.15))
        k_s = k_yr/(365*24*3600)
        return (100*self.host.E_synth*k_s)/cutoff_pc

    """ Tijhuis calcualtion [QE Document] """

    def getTijuisME(self):
        drym_ME = 4.5*math.exp((self.host.locale.env.T**(-1)-298**(-1))*((-6.94*10000)/8.31))
        return (1000/3600)*(self.host.dry_mass/0.026)*drym_ME
