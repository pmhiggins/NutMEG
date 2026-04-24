Creating an ecosystem
=====================

NutMEG is designed for solving problems involving multiple organisms. To simulate this, we collect them together in a shared Reactor in an Ecosystem. First, a community of similar organisms can be defined in an :any:`OrganismPopulation` then any number of OrganismPopulations can be initialised inside an :any:`Ecosystem`. Continuing with the Reactor and BaseOrganism's created on previous pages:

..code::

  import NutMEG.core as nmc


  MGPop = nmc.OrganismPoplation(base_MG, 100) # 100 methanogens

  # You can equivalently use nmc.OrgPop
  SRPop = nmc.OrgPop(base_SR, 1000) # 1000 sulfate reducers

  ES = nmc.Ecosystem(R, [MGPop, SRPop]) # R is the reactor we built earlier

The EcoSystem is used as the 'hub' where we run `NutMEG.models` to solve biogeochemical and astrobiological problems. Continue reading to check out some simple examples of what can be explored!
