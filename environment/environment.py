"""

This is the environment module, kept seperate from the reaction module for clarity,
though the two really do depend on each other. For now, any changes that a reaction
makes on its environment will be contained in the reaction module and simply update
our environment.

@author P M Higgins
@date 10/2017
@version 0.0.1
"""

class environment:
    """
    Class to keep our environment parameters all in one place for use
    in various reactions. In its current form, this is a class of
    convenience i.e. it doesn't really do anything --- just stores all
    the physical variables of the world enclosed by volume V takes
    place in. This may well change in the future though.

    SI units are used throughout.
    """

    T = None # average temperature of container/vessel in K
    V = None # total volume of our container/vessel in m^3
    P = None # total pressure of container/vessel in N/m^2


    # RTP conditions
    T_RTP = 298.15 # standard temperature
    V_RTP = None # total volue of vessel at RTP
    P_RTP = None # total pressure of vessel at RTP (usually 1 bar,
      # but will vary for aqueous environments.

    def __init__(self, T=298.15, V=0.001, P=101325.0, V_RTP=None, P_RTP=None):
        # if nothing is passed, we assume RTP temperature,
          # atmospheric pressure, and
          # a volume of 1 L (10^-3 m^-3)
        self.T=float(T)
        self.V=float(V)
        self.P=float(P)
        self.V_RTP=V_RTP
        self.P_RTP=P_RTP
