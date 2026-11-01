Bioenergetics 101
=================

NutMEG is built on the idea of energetic habitability. Energy is the universal
currency for life as we understand it, and while metabolic processes vary, they
will have a net energetic cost or yield which is quantifiable via standard
chemical thermodynamics. By identifying energy sources and processes
corresponding to maintenance (e.g. survival in adverse conditions) and growth
we can assess the habitability of various environments. [1]_

If the energetic availability in an environment outweighs the energetic cost of
surviving there, we say it is energetically habitable. Where possible, NutMEG
estimates microbial growth via an approach where the
energetic input from metabolism (:math:`P_{S}`) is directed toward microbial maintenance
(:math:`P_{M}`), and growth (:math:`P_{G}`) - hence the acronym from *Nutrients,
Maintenance, Energy and Growth*.

.. math::
   P_{G} = P_{S} - P_{M}

The maintenance power :math:`P_{M}` reflects the total cost of microbial maintenance
processes - ones which are necessary for survival but do not directly contribute
to growth. These could include maintaining a specific internal pH
(:math:`P_{pH}`), repairing biomacromolecules as they break down with
temperature (:math:`P_{T}`), and defending against adverse salinity
(:math:`P_{SAL}`) to give a few examples. Mathematically it takes this form:

.. math::
   P_{main} = P_{pH} + P_{T} + P_{SAL} + ...

The nutrient uptake efficiency is more complex to compute. It represents the effect
of limited availability of carbon, hydrogen, nitrogen, oxygen, phosphorus, or
sulfur (CHNOPS) elements on the amount of biomass an organism can actually make
per unit time. [2]_


-----

.. [1] Hoehler (2007) *Astrobiology* doi: `10.1089/ast.2006.0095 <hhttps://doi.org/10.1089/ast.2006.0095>`_
.. [2] Higgins & Cockell (2020), *J. R. Soc. Interface* doi: `10.1098/rsif.2020.0588 <hhttps://doi.org/10.1098/rsif.2020.0588>`_,
