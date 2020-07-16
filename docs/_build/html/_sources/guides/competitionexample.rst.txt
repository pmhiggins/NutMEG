Simple example: simulating competition
======================================

An example of a simple simulation can be found in the `Competition_Example <https://github.com/pmhiggins/NutMEG-Implementations/Competition_Example>`_ directory of `NutMEG-Implementations <https://github.com/pmhiggins/NutMEG-Implementations>`_, in which we explore how different environmental and microbial parameters can affect competition between some (entirely fictional) methanogens and sulfate reducers. In this guide we'll go in-depth into how this code works for users who may not be used to object-oriented programming in Python.

``setup_methanogenesis()`` and ``setup_sulfatereduction()``
-----------------------------------------------------------
These methods generate ``reaction`` objects to represent the two metabolisms we're looking at. Each method initialises several ``reagent`` objects to be included in the overall reaction. Note how the ``reactor`` is passed here, this is to ensure that all the ``reagent`` s and ``reaction`` s point to the same ``environment`` object. You could also simply pass the ``environment`` object.

Note how the rate constant at RTP (room temperature and pressure) is passed as an optional argument. This is one parameter (or alternatively, the rate constant at the current environment, ``rate_constant_env``) which is necessary for the kinetic calculations (see `respirator <../source/NutMEG.culture.base_organism.html>`_). You could not pass a value, and the default k_RTP will be that for an optimally growing methanogen at 300 K. 

``initial_conditions()``
------------------------
This method sets up the ``reactor`` by populating it with ``reagent`` s. You may notice some of these are the same reagents defined in the two methods discussed above but as they're being initialised in separate functions they are completely different objects. You may also notice that these reagents are being initialised with more information than the others too, including their activities and concentrations. That;s becaus these are the key ones we'll be using in the simulation. We'll see how everything is unified in the next section.

Some of the reagents added have nothing to do with the metabolisms in ``setup_methanogenesis()`` and ``setup_sulfatereduction()`` like NH_3. This is because NutMEG also monitors nutrient availability when it predicts microbial behaviour. The majority of biomass is made up of so-called CHNOPS elements, and nothing will grow in NutMEG without some of each (unless you built your own organism with a different biomolecular structure to the defaults in NutMEG).

``simulate_competition()``
--------------------------
This method sets up and runs a growth prediction simulation. First, a reactor is created:

.. code::

    R = nm.reactor('reactor1', workoutID=False, pH=7.0, dbpath=dbpath, **reactor_changes)

here ``'reactor1'`` is its name, which tells the database to have a specific table for this reactor type. Passing ``workoutID=False`` means it doesn't save parameters to the database yet. We'll do so once its fully populated. ``**reactorchanges`` passes the dictionary ``reactorchanges`` as keyword arguments to the ``reactor.__init__`` method. Next, we set up the reactor's composition using ``initial_conditions``. Then we set up the organisms, as ``horde`` objects:

.. code::

    H = nm.horde('Methanogen', R, setup_methanogenesis(R, k_RTP=k[0]), 10,
      Tdef='None', dbpath=dbpath, **methanogen_changes)

    H2 = nm.horde('SulfateReducer', R, setup_sulfatereduction(R, k_RTP=k[1]), 10,
      Tdef='None', dbpath=dbpath, **SR_changes)

They have separate names, so they'll also have their own tables in the database. They are both put into the same reactor ``R`` and have a unique metabolic reaction. When they're added to ``R``, the reagents are unified to those already populating ``R`` to avoid accidentally 'doubling-up' your reagents. ``Tdef`` shows which adaptations against temperature to use, here we'll ignore it and put in our own maintenance powers.

Now that ``R`` is fully populated, we can give it an entry in the database:

.. code::

    R.dbh.workoutID()

``ecosystem`` requires organisms and hordes to be in ``culture`` objects for easier analysis of them. Then, the culture and reactor just need to be put in an ecosystem together, and a growth prediction can be performed with ``predict_growth``:

.. code::

    Cu = nm.culture(hordes=[H2, H])
    ES = nm.ecosystem(R, Cu, dbpath=dbpath)
    ES.predict_growth(tmax=2e8)


``orgcurves()`` and ``compcurves()``
-------------------------------------
These methods extract the simulation data from the database and plot them using the ``NutMEG.plotter`` module which specifically uses matplotlib as a basis, but you could easily use another python plotting package.

Using a SimID returned by ``simulate_competition()``, ``orgcurve()`` extracts and plots the parameter ``'param'`` passed from the database with time. For a list of possible values of ``'param'`` check out the guide on database structure.

Similarly, ``compcurve()`` plots the concentration of select reagents with time.
