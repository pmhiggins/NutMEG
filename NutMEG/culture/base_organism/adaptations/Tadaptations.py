"""

Class for calculating adaptations to adverse T conditions

"""
import math

class Tadaptations:

    TOM_1 = [ 3.75274464e-10, -6.28614763e-07,  4.20934430e-04,
      -1.40873711e-01, 2.36260582e+01, -1.60782148e+03]
    TOM_05 = [ 5.47545848e-10, -9.36824600e-07,  6.41040964e-04,
      -2.19319054e-01, 3.75813051e+01, -2.59970377e+03]
    TOM_15 = [ 1.40807504e-10, -2.37559221e-07,  1.60324386e-04,
      -5.41193986e-02, 9.19861608e+00, -6.48600497e+02]


    def __init__(self, host):
        self.host = host # the host oragnism and its assosiated environment.

    # Lever calculation [QE Document]
    def getLeverME(self, cutoff_pc=10):
        """Power cost due to temperature according to Lever et al 2015
        @param cutoff_pc the % racemization of AA in proteins at which it is
        replaced."""

        k_yr = 0.00012*math.exp(0.10174*(self.host.locale.env.T-273.15))
        k_s = k_yr/(365*24*3600)
        return (100*self.host.E_synth*k_s)/cutoff_pc

    # Tijhuis calcualtion [QE Document]
    def getTijhuisME(self):
        """Power cost due to temperature according to Tijhuis 1993"""
        drym_ME = 4.5*math.exp((self.host.locale.env.T**(-1)-298**(-1))*((-6.94*10000)/8.31))
        return (1000/3600)*(self.host.dry_mass/0.026)*drym_ME

    def getTijhuisAerobe(self):
        """Power cost due to temperature according to Tijhuis 1993"""
        drym_ME = 5.7*math.exp((self.host.locale.env.T**(-1)-298**(-1))*((-6.94*10000)/8.31))
        return (1000/3600)*(self.host.dry_mass/0.026)*drym_ME

    def getTijhuisAnaerobe(self):
        """Power cost due to temperature according to Tijhuis 1993"""
        drym_ME = 3.3*math.exp((self.host.locale.env.T**(-1)-298**(-1))*((-6.94*10000)/8.31))
        return (1000/3600)*(self.host.dry_mass/0.026)*drym_ME

    def getTOM(self):
        """Power cost due to temperatre for methanogens from Higgins & Cockell 2020"""
        HCpolyfit = None
        if self.host.respiration.n_ATP >= 1.5:
            HCpolyfit = self.TOM_15
        elif self.host.respiration.n_ATP <= 0.5:
            HCpolyfit = self.TOM_05
        else:
            HCpolyfit = self.TOM_1

        MP = 10**sum([x*(self.host.locale.env.T**(5-i)) for i, x in enumerate(HCpolyfit)])
        return MP * self.host.dry_mass / (300*3.4412868852915668e-18 )
