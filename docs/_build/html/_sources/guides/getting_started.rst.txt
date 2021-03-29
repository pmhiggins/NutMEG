
Getting Started
================

This short guide will help you get NutMEG up and running. You can download the
full package from `github <https:github.com/pmhiggins/NutMEG>`_.

Installing Dependencies
-----------------------
NutMEG has been built and is running on python 3.7. Previous versions of python
may be stable but have not been tested. Users have noted errors occurring 
for systems running on python 3.8 and above. We are aware of this and will be 
pushing an unpdate in the near future.

In order to perform certain chemical calculations, NutMEG has a dependency on
the `reaktoro <http://en.wikipedia.org/wiki/Hyperlink>`_ package. ``reaktoro``
is hosted on conda-forge, so we recommend using NutMEG in a conda environment
for ease of use.

By far the simplest route to get up and running would be to install
`miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ then prepare to
use NutMEG in a ``conda`` environment. In a terminal, write:

.. code::

    conda config --append channels conda-forge
    conda install reaktoro

NutMEG also depends on a few typical python libraries including ``numpy``,
``pandas``, and ``matplotlib``. If you don't have these installed in your
conda environment navigate to your NutMEG directory and run:

.. code::

    pip install -r requirements.txt


Structuring NutMEG Projects
---------------------------
While NutMEG is in development and not widely released, it will simply have to
sit in project folders or somewhere you can append to ``sys.path``. For example,
you could structure your project folder like so:

.. code::

    your_project/
    |-- NutMEG/
    `-- your_project_code.py
