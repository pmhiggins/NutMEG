Creating a reactor
==================

The :any:`Reactor` class represents the physicochemical environment of a NutMEG simulation. Reactors should be initialised first, then can be populated with
:any:`Reagent` or :any:`Reaction` instances you could want to capture. For simulating organisms, the local Reactor must host an overall metabolic equation (more on that on the next page).

Creating reactors from scratch
------------------------------

Below is a minimal example setting up a :any:`Reactor` that can model H2 methanogenesis and sulfate reduction.

.. code:: python

  import NutMEG.core as nmc

  R = nmc.Reactor(pH=7) # default is 1 L at RTP, so T=298 K, P=1e5 Pa, 1L H2O

  # introduce some reagents
  # The string identifier should point to a species present in Reactor's
  # thermodb attribute, which is SUPCRT07 by default.

  CO2 = nmc.Reagent('CO2(aq)', R, phase='aq', amount=(0.001, 'mol'))
  H2aq = nmc.Reagent('H2(aq)', R, phase='aq', amount=(0.001, 'mol'))
  CH4aq = nmc.Reagent('Methane(aq)', R, phase='g', amount=(1e-7, 'mol'))
  SO4 = nmc.Reagent('SO4-2', R, phase='aq', charge=-2, amount=(0.001, 'mol') )
  HS = nmc.Reagent('HS-', R, phase='aq', charge=-1, amount=(1e-7, 'mol'))
  H2O = nmc.Reagent('H2O(aq)', R, phase='l', activity=1.)

  # add the methanogenesis reaction
  thermalMG = nmc.Reaction(R, {'CO2(aq)':1, 'H2(aq)':4}, {'Methane(aq)':1, 'H2O(aq)':2})

  # add sulfate reduction
  thermalSR = nmc.Reaction(R, {'H2(aq)':4, 'SO4-2':1, 'H+':1}, {'HS-':1, 'H2O(aq)':4})

Reactors can perform rudimentary chemical reactions with the ``perform_reaction`` function:

.. code::

  R.perform_reaction(thermalMG.equation, 1e-5) # react 1e-5 moles of methanogenesis

They can also be simulated as open systems by assigning the ``composition_inputs`` attribute:

.. code::

  # every second, 1mmol of CO2 is added and 1mm of CH4 is removed
  R.composition_inputs = {'CO2': 1e-3, 'Methane(aq)': -1e-3}

This and reactions are dynamically updated by higher-level time-integration modules in NutMEG we will explore in the examples.

Importing from reaktoro
-----------------------

If you use reaktoro for chemical modelling, a NutMEG reactor can easily be built from a reaktoro state! The syntax is as below

.. code::

  # imports all aqueous species and their thermodynamic properties
  # from a reaktoro chemical system.
  R1 = nmc.Reactor.from_rto_state(your_ChemicalState, your_ChemicalProps)

  # if you iterate between reaktoro and NutMEG calculations, it is also
  # possible to update an existing reactor from a reaktoro state:
  R1.update_from_rto_state(your_ChemicalState, your_ChemicalProps)
