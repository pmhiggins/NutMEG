
Installation
============

Installing with conda and pip
-----------------------------

NutMEG has been built and is running on python 3.12. Previous versions of python
may be stable but have not been tested.

In order to perform certain chemical calculations, NutMEG has a dependency on
the `reaktoro <http://reaktoro.org/>`_ package. ``reaktoro``
is hosted on conda-forge, so we recommend using NutMEG in a conda environment
for ease of use.

To get up and running, install
`miniconda <https://docs.conda.io/en/latest/miniconda.html>`_
(or your favorite conda distributor) then prepare to use NutMEG in a ``conda``
environment. The following termianl commands create a conda environment called
``nm`` (use any name you wish), activates it, and installs NutMEG there:

.. code::

    conda config --add channels conda-forge
    conda create --name nm reaktoro
    conda activate nm
    pip install NutMEG-astrobiology

You can now ``import NutMEG`` provided your conda environment is activated.

For minor bugfixes while v2 is in beta (but not major ones, causing a new release),
NutMEG can also be  be installed/updated from  `github <https:github.com/pmhiggins/NutMEG/tree/nm_v2>`_. via:

..code::

    pip install git+https://github.com/pmhiggins/NutMEG.git@nm_v2

Installing via `pip` should also install NutMEG's additional required packages.
If not, or there is an error, the requirements are listed below:

.. code::
    "numpy>=1.26.4",
    "pandas>=2.2.3",
    "uncertainties>=3.2.2",
    "PyYAML",
    "reaktoro>=2.12.5",

.. note ::

    NutMEG v1 only worked with reaktoro v1.x, and NutMEG v2 only works with
    reaktoro v2.x. Make sure you install the correct version.
