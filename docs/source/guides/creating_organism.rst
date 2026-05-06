Creating an organism
====================

The :any:`BaseOrganism` class represents a model organism and collects its organismic properties. Properly setting up a BaseOrganism is more involved than it was in NutMEG v1 to allow for maximal flexibilty in what assumptions we want to make about organisms.

A minimal initialisation requires the overall metabolic reaction the organism performs, and it's maximum metabolic rate and growth rate and could take the form:

.. code::
  base_MG = nmc.BaseOrganism('Methanogen', thermalMG)
  base_MG.metabolism.set_max_rate(1e-17) # mol/s/cell
  base_MG.growth.set_max_rate(float('inf')) #/s

This would allow the organism to grow as fast as it can practically access energy via its metabolic reaction. :any:`BaseOrganism` computes its bioenergetic and biokinetic parameters using three helper classes: :any:`Metaboliser`, :any:`Maintainer`, and :any:`Grower`. A more involved organism initialisation can directly create this attributes for more control. Let's make a sulfate reducer to compare:

.. code::

  import NutMEG.models as nmm

  SR_MetRate = nmm.base_rate_models.FirstOrderChemical(5e-4)
  SR_GroRate = nmm.base_rate_models.ConstantRate(float('inf'))

  # we can also define properties that affect the metabolic rate or growth rate
  SR_Monod = nmm.forcing_factors.Monod('H2(aq)', 1e-3)
  SR_BE = nmm.forcing_factors.Bioenergetic(n_ATP=1.)
  SR_Phosphate = nmm.forcing_factors.Monod('PO4', 1e-3)

  # now build our model sulfate reducer
  base_SR = nmc.BaseOrganism(
    'SulfateReducer',
    nmc.org.Metaboliser(thermalSR, base_rate=SR_MetRate, forcing_factors=[SR_Monod, SR_BE]),
    growth = nmc.org.Grower(base_rate=SR_GroRate, forcing_factors=[SR_Phosphate])
  )

Here, the sulfate reducers rate can be controlled by both the available [H2] (via a Monod expression), Bioenergetics and its ATP energy yield, and the available [PO4] (via a seperate Monod expression). Check out the :any:``NutMEG.models.organism`` subpackage for a variety of rate-limiting models to try!
