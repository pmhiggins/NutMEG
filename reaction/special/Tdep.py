"""
This is the Tdep submodule of the reaction module. It emphasises the
temperature dependence of reactions. For now we have only the Tdep class,
which will correct reaction for discrepancies in H(T) and S(T), if known.
We may add further classes later for higher temperature fidelity.

Usually, a standard reaction with reaktoro will do. However, if we don't
have much data this could come in handy.

@author P M Higgins
@version 0.1
"""

import sys
sys.path.append('..')
from NutMEG.reaction import reaction
import numpy as np
import math

class Tdep(reaction):
    """Class for temperature dependent free energy calculations.

    For use in a scenario where we do not have anything available in a
    thermochemical database but do have the RTP values and specific heat
    capacities of the reagents. Then, there are a few routes to \DeltaG(T).
    """


    """
        INITIALISATION
    """

    # initialise as a reaction, but only wih the relevant attributes.
    # add other things in as we go along.
    def __init__(self, reactants, products, env, equilibrium=False,
      *args, **kwargs):
            reaction.reaction.__init__(self, reactants, products, env,
              equilibrium=equilibrium, *args, **kwargs)


    """
        SCALING REACTION PARAMETERS FROM RTP
    """


    # Looks like the integral is working but you might want to test more thoroughly

    def update_std_molar_enthalpy_of_reaction_T(self):
        """Update the standard molar enthalpy at this temperature.

        Use reactant paramters to get the non-RTP \Delta H(T)^o,
        providing Cp takes integer polynomial form.
        This calculation is valid for ideal gases and infinitely
        diluted solutions.

        """
        # Calculate the enthaply of formation for each reagent at the required temperature
        HoR =0.
        for p, mr in self.products.items():
            if p.Cp_T_poly is None:
                raise ValueError("You have not assigned specific heat "
                  "capacities to all of your products!")
            integral = np.polyint(p.Cp_T_poly)
              # ^ no need for constants as we're going to add in the limits
            HoR = HoR + (p.std_formation_enthalpy_RTP
              + (mr * (integral(self.env.T) - integral(self.env.T_RTP))))
        for r, mr in self.reactants.items():
            if r.Cp_T_poly is None:
                raise ValueError("You have not assigned specific heat "
                  "capacities to all of your reactants!")
            integral = np.polyint(r.Cp_T_poly)
            HoR = HoR - (r.std_formation_enthalpy_RTP
              + (mr * (integral(self.env.T) - integral(self.env.T_RTP))))
        self.std_molar_enthalpy = HoR


    def update_std_molar_entropy_of_reaction_T(self):
        """Update the standard entropy of at this temperature.

        Use reactant paramters to get the non-RTP \Delta S(T)^o,
        providing Cp takes integer polynomial form.
        This calculation is valid for ideal gases and infinitely
        diluted solutions.

        """

        SoR = 0.
        for p, mr in self.products.items():

            Cp_TbyTarr = [] # to calculate the Cp/T polynomial

            if p.Cp_T_poly is None:
                raise ValueError("You have not assigned specific heat "
                "capacities to all of your products!")
            if p.Cp_T_poly == np.poly1d([0]):
                SoR += p.std_formation_entropy_298
                # if we've inititalised the heat capacity as zero,
                  # skip this reagent
                continue

            i=0
            invTterm = 0.
            for cp in p.Cp_T.c:
                i = i+1
                if len(p.Cp_T_poly.c) != i:
                    # ignore the last term. With the 1/T any constants will
                      # cancel when we apply the limits anyway.
                    Cp_TbyTarr.append(cp)
                    # ^ by the end of the loop this is Cp/T in polynomail form.
            integral = np.polyint(np.poly1d(Cp_TbyTarr))
            SoR = SoR + (p.std_formation_entropy_RTP
              + (mr * (integral(self.env.T) - integral(self.env.T_RTP)
              + math.log(self.env.T/self.env.T_RTP))))
        for r in self.reactants:
            Cp_TbyTarr = []
            if r.Cp_T_poly is None:
                raise ValueError("You have not assigned specific heat "
                  "capacities to all of your reactants!")
            if p.Cp_T_poly == np.poly1d([0]):
                SoR -= r.std_formation_entropy_RTP
                continue
                # if we've inititalised the heat capacity as zero,
                  # skip this reagent
            i=0
            invTterm = 0.
            for cp in r.Cp_T_poly.c:
                i = i+1
                if len(r.Cp_T_poly.c) != i:
                    # ignore the last term. With the 1/T any constants will
                      #cancel when we apply the limits anyway.
                    Cp_TbyTarr.append(cp)
                    # ^ by the end of the loop this is Cp/T in polynomail form.
            integral = np.polyint(np.poly1d(Cp_TbyTarr))
            SoR = SoR - (r.std_formation_entropy_RTP
              + (mr * (integral(self.env.T) - integral(self.env.T_RTP)
              + math.log(self.env.T/self.env.T_RTP))))
        self.std_molar_entropy = SoR




    # Then, K_P. However, to realise this we need H as a function of T.
    # We shouldn't need it yet so I'll just leave it for now.


    # A way on integrating other forms of Cp would be useful. This can
      # be done through scipy.integrate.quad. The only problem is that
      # we'd have to explicitly define each Cp equation somewhere in the code.

    # If I build a database for it to use in the future, then they can
    # be saved in there somewhere...
