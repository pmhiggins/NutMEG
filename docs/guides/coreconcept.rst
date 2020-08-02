The core concept
================

NutMEG is built on the idea of energetic habitability. Energy is the universal
currency for life as we understand it. While metabolic processes vary, they
will have a net energetic cost or yield which is quantifiable via standard
chemical thermodynamics. By identifying energy sources and processes
corresponding to maintenance (e.g. survival in adverse conditions) and growth
we can assess the habitability of various environments.

If the energetic availability in an environment outweighs the energetic cost of
surviving there, we say it is energetically habitable. NutMEG estimates
microbial growth via an approach where efficiencies are applied to the
energetic input from metabolism (:math:`P_{S}`) corresponding to  microbial maintenance
(:math:`\epsilon_{M}`), and nutrient uptake (:math:`\epsilon_{UT}`). Any
leftover energy can be directed into biomass synthesis (:math:`P_{G}`).

.. math::

    P_{G} = \epsilon_{UT} \epsilon_{M} P_{S}

The maintenance efficiency reflects the total cost of microbial maintenance
processes - ones which are necessary for survival but do not directly contribute
to growth. These could include maintaining a specific internal pH
(:math:`P_{pH}`), repairing biomacromolecules as they break down with
temperature (:math:`P_{T}`), and defending against adverse salinity
(:math:`P_{SAL}`) to give a few examples. Mathematically it takes this form:
