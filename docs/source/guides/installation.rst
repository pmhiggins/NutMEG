
Installation
============

Installing with pip
-------------------

...is coming soon!

Until v2 releases, NutMEG must be built from  `github <https:github.com/pmhiggins/NutMEG/tree/nm_v2>`_.


Installing Dependencies
-----------------------
NutMEG has been built and is running on python 3.7-12. Previous versions of python
may be stable but have not been tested.

In order to perform certain chemical calculations, NutMEG has a dependency on
the `reaktoro <http://reaktoro.org/>`_ package. ``reaktoro``
is hosted on conda-forge, so we recommend using NutMEG in a conda environment
for ease of use.

By far the simplest route to get up and running would be to install
`miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ then prepare to
use NutMEG in a ``conda`` environment. In a terminal, write:

.. code::

    conda config --append channels conda-forge
    conda install reaktoro

.. note ::

    NutMEG v1 only worked with reaktoro v1.x, and NutMEG v2 only works with
    reaktoro v2.x. Make sure you install the correct version.

NutMEG also depends on a few typical python libraries including ``numpy``,
``pandas``, and ``matplotlib``. If you don't have these installed in your
conda environment navigate to your NutMEG directory and run:

.. code::

    pip install -r requirements.txt


Structuring NutMEG Projects
---------------------------
Until NutMEG has a PyPI release, there are two options to use it from source.

1. Place the repository in project folders or somewhere you can append to ``sys.path``. For example,
you could structure your project folder like so:

.. code::

    your_project/
    |-- NutMEG/
    `-- your_project_code.py

2. Permanently add it to your ``sys.path`` so you can more
easily import NutMEG anywhere. Navigate to the ``site-packages`` folder of
your conda environment and add a file ``usercustomize.py`` containing:

..code::

  import sys
  sys.path.extend(['/path/to/NutMEG'])
