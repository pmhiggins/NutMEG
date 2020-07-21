import NutMEG

class organism(NutMEG.base_organism):
    """Class for an individual organism, which grows and splits
    independently.
    """


    def __init__(self, name, locale, metabolism,
      maintenance=None,
      CHNOPS=None,
      mass=1.e-16,
      E_synth=None,
      volume=1.e-18):
        NutMEG.base_organism.__init__(self, name=name, locale=locale,
          metabolism=metabolism,
          maintenance=maintenance,
          CHNOPS=CHNOPS,
          mass=mass,
          E_synth=E_synth,
          volume=volume)



    def take_step_lite(self, t):
        """Increment the organism's life by time t.

        Updates the organism's parameters based on its mortality,
        environment, metabolism, etc.
        """
        if self.isactive:
            self.age += t
            # EVERYTHING COMMENTED BELOW HERE HAS MOVED INTO COLONY LITE

            #power supplies updated in colony_lite

            # update energy used
            self.E_store += self.maintenance.get_P_store()*t
            # E_growth and CHNOPS done in colony_lite

            # update volume as it grows
            self.volume = self.base_volume*(1.0 + (self.E_growth/self.E_synth))

            #reaction performed in colony_lite

            if self.E_growth > self.E_synth:
                # split the cell
                if round(self.E_growth/self.E_synth)>1:
                    raise ValueError('Your timestep is too long, '
                      + self.name + ' is splitting multiple times per '
                      + 'timestep!')
                else:
                    self.issplitting = True
                    # the organism will be split by the colony

        if self.age > self.life_span:
            # kill the organism. It remains as biomass,
            # but cannot divide or metabolise
            self.isactive=False
