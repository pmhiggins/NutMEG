Core framework
==============

NutMEG is arranged into three subpackages:

* :any:`NutMEG.core` contains the building-block objects for creating a wide variety of computational scenarios. They are designed to be as simple and flexible as possible, and can be extended for specific applications and presets.

* :any:`NutMEG.models` contains classes that for models used by :any:`NutMEG.core` objects. For example, a :any:`BaseOrganism` requires a growth_model describing its growth behaviour. Classes that can faciliate that are found in the  :any:`NutMEG.models.organism.growth_models` subpackage, like :any:`Monod`, or :any:`BioenergeticGrowthModel`. They can also be user-defined.

* :any:`NutMEG.presets` provides builders for a specific type of organism or environment. For example, :any:`Enceladus` [1]_ creates an Enceladus-like chemical environment, and :any:`TypicalOptimalMethanogen` creates methanogenic organism population as defined in ref [2]_.



The basics of ``NutMEG.core``
-----------------------------


``NutMEG.core.ecosystem``
^^^^^^^^^^^^^^^^^^^^^^^^^

Most presets and applications interface with an :any:`Ecosystem`. It contains the local chemical environment which is summarised by a :any:`Reactor` object, and any number of :any:`Resident` objects that can interact with that environment through time.


The ``NutMEG.core.reactor`` submodule
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* :any:`Reactor` objects monitor the physicochemical system, such as it's composition, size, temperature, and pressure. Under default settings, :any:`Reactor` is valid for well-mixed near-equilibrium or kinetically quiescent conditions where the composition does not appreciably change unless we introduce something to disrupt the system such as an :any:`OrganismPopulation`. Simple environmental disequilibria can be introduced, such as a constant fresh supply of a :any:`Reagent` or simple modelled chemical :any:`Reaction`.

* :any:`Reagent` objects store and updates properties of chemical species such as thermodynamic parameters, activities, molar quantities. Thermodynamic data at a range of temperatures and pressures can be calculated from chemical databases.

* :any:`Reaction` objects model chemical reactions by interacting with ``Reagents``. Both biotic and abiotic chemical interactions with the :any:`Reactor` are possible.

The ``NutMEG.core.base_organism`` submodule
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This submodule is for building model organisms.

* :any:`BaseOrganism`

* :any:`Metaboliser`

* :any:`Maintainer`

* :any:`Grower`


The ``NutMEG.core.resident`` submodule
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* :any:`OrganismPopulation`



---

.. [1] Higgins et al (2021) *JGR:Planets*; doi: `10.1029/2021JE006951 <https://doi.org/10.1029/2021JE006951>`_,
.. [2]  Higgins & Cockell (2020), *J. R. Soc. Interface* doi: `10.1098/rsif.2020.0588 <hhttps://doi.org/10.1098/rsif.2020.0588>`_,
