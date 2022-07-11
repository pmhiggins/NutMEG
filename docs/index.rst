.. NutMEG documentation master file, created by
   sphinx-quickstart on Wed May 13 15:18:38 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


NutMEG
======

*Nutrients, Maintenance, Energy and Growth*

Welcome to NutMEG's documentation! `NutMEG <https://github.com/pmhiggins/NutMEG>`_ is a python module for predicting the
growth behaviour of microbial organisms in astrobiology. It's designed to
estimate the energetic availability through local chemistry and the cost of
defending against extreme or adverse conditions to work out whether an
environment could be habitable and how biology could behave there.

NutMEG is still a work-in-progress, but as of its most recent release it can be
used to create tentative growth predictions of known organisms in known environments.
Check out the `competition example <guides/competitionexample.html>`_ for an idea of
what it can do at this stage.

Naturally, the above means this documentation is also a work-in-progress. In the meantime
the most thorough summary of NutMEG and its capapbilities can be found in `this thesis <http://dx.doi.org/10.7488/era/2078>`_,
specifically Chapters 3 (for theory), 4 (for code design), and the appendix for bonus code examples.


.. toctree::
    :maxdepth: 1
    :caption: Guides

    guides/coreconcept
    guides/getting_started
    guides/class_structure
    guides/database_structure
    guides/growthalgorithm
    guides/implementations
    guides/competitionexample
    guides/citing_nutmeg

.. toctree::
  :maxdepth: 1
  :caption: contact

  contact/get_in_touch

.. toctree::
    :maxdepth: 3
    :caption: API

    source/modules



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
