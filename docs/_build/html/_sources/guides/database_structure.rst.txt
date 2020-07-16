Database Structure
==================

Overview
--------
NutMEG simulations can become expensive, particularly if you're performing calculations
over a wide range of parameters or working over long timescales. To prevent you having to rerun
simulations NutMEG saves a copy of the state of the ecosystem every 100 steps up until the final
step by default to a named SQL database. This is by default named 'NutMEG_db' and placed in the
working directory, but can be changed by specifying ``dbpath``, a common keyword argument for objects in
NutMEG, or changing ``std_dbpath`` in ``NutMEG.util.NutMEGparams.py``.

Three key identifiers are used throughout NutMEG's database: SimID, OrgID, and LocID.
SimID is the primary key of the Summary table, and each one links to a table with the
name ``'FullResults_Sim_'+SimID``.

LocID identifies ``reactor``-like objects and is a primary key of the Reactor table.
Similarly OrgID identifies ``base_organism``-like objects and is a primary key of
the organism table. Of course, organisms and reactors can take different forms and
have different attributes depending on the level of complexity you want to achieve
with your simulations and whether you're using ``saved_systems`` or saved_organisms.
Instance variables of reactors and organisms are store in the NutMEG database so as to
tell if a simulation has ben performed before or not.
Each instance of a ``reactor``-like or ``base_organism``-like object has a ``name``
attribute. This variable tells NutMEG which table the objects parameters should be
saved to, handy if you've included your own ``output`` helper class. It also assigns
the prefix of the relevant OrgID or LocID.


Tables
-------

Below is the structure of some of the tables, populated with some example organisms, reactors, and simulations.

**Summary Table: static results for each simulation**


+-----------+---------------------------+---------------------+----------------+--------------+-----+
| SimID     | OrgIDs                    | LocID               | FinBM          | PeakGR       | etc |
+===========+===========================+=====================+================+==============+=====+
| 1_010120  | (’Methanogen_1_010120’,)  | Enceladus_1_010120  | :math:`10^{6}` | :math:`10`   |     |
+-----------+---------------------------+---------------------+----------------+--------------+-----+
| 2_010120  | (’Methanogen_2_010120’,)  | Enceladus_1_010120  | :math:`10^{5}` | :math:`7`    |     |
+-----------+---------------------------+---------------------+----------------+--------------+-----+


 **FullResults_Sim_1_010120 Table:**


+------+-----------------+------------------+-----------------+------+
| time | composition     | no_alive_Methan\ | EnergyAvai\     | etc. |
|      |                 | ogen_1_010120    | lable_Methan... |      |
+======+=================+==================+=================+======+
| 100  | CO2:0.01,       | 550              | 141000          |      |
|      | CH4:0.001...    |                  |                 |      |
+------+-----------------+------------------+-----------------+------+
| 200  | CO2:0.009,      | 600              | 140900          |      |
|      | CH4:0.0011...   |                  |                 |      |
+------+-----------------+------------------+-----------------+------+


**Organism table. Match up organisms to metabolisms so we
know which table to pull up for them**


===================== ===============
OrgID                 Type
===================== ===============
Methanogen_1_010120   Methanogen
SulfateR_1_010120     Sulfate Reducer
Methanogen_2_010120   Methanogen
Methanogen_1_020120   Methanogen
===================== ===============

**Methanogen Table:**


+----------+----------+----------+----------+----------+------+
| OrgID    | Res\     | Esynth   | dry mass | Tdef     | etc. |
|          | piration |          |          |          |      |
+==========+==========+==========+==========+==========+======+
| Metha\   | CO2 + \  | :math:`\ | :math:`\ | Lever10pc|      |
| nogen_1\ | 4H2 -> \ | 7.9\     | 1\times\ |          |      |
| _020120  | CH4 + \  | \times\  | 10^{-13}`|          |      |
|          | 2H2O     | 10^{-11}`|          |          |      |
+----------+----------+----------+----------+----------+------+

 **Reactor table. Match up locales to named systems so we
 know which table to pull up for them**

==================== =========
LocID                Env Type
==================== =========
Enceladus_1_010120   Enceladus
Mars_1_010120        Mars
Enceladus_2_020120   Enceladus
aqueous_24_030220    Aqueous
==================== =========



**Enceladus Table:**

+----------------------+-------+----------+-------------+------------------------+------+
| LocID                | Depth | Pressure | Temperature | Composition            | etc. |
+======================+=======+==========+=============+========================+======+
| Enceladus_1_010120   | 8.5   | 40       | 290         | CO2:0.01, CH4:0.001... |      |
+----------------------+-------+----------+-------------+------------------------+------+


Extracting data for analysis
----------------------------

.. note :: A full guide on this is forthcoming...

Generally you can extract any simulation parameters using the static methods found in
`ecosystem_dbhelper <../source/NutMEG.html#module-NutMEG.ecosystem_dbhelper>`_
. The method ``extract_param_db_Sim()`` will extract parameters for the given
SimID pertaining to the given orgID, if specified, and return them. These could be
entries in Summary (e.g. PeakGR for the maximum growth rate in the simulation), or entries
in the FullResults_Sim\_ table (e.g. no_alive for the active current population throughout
the simulation).

You can't extract specific parameters from the organism and reactor relevant tables,
but you can initialise an identical object using the ``r_from_db()`` or ``bo_from_db()`` class methods --- think of this as similar to pickling.
