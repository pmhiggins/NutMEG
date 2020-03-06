import sys
#sys.path.append("../lib/rto/lib/python2.7/site-packages") # local import of reaktoro
import os
from reaktoro import *
import numpy as np
import math

this_dir = os.path.dirname(__file__)


class BioMolecule:
	"""

	class for storing amino acid properties, or any other useful biological molecule really

	"""
	name = ''
	Mr = 0. # molecular mass
	conc_mol_per_l = None	# mol/l
	conc_mol_per_cell = None # mol/cell

	gamma = 1.		# activity coefficient

	std_formation_gibbs = None		# J/mol
	std_formation_R = None #J/mol

	#booleans of state
	thermo = True # whether we have the thermodynamic data available

	copies_per_cell = None
	frequency = None
	probrange = None

	#generic counter to use how you please
	counter = 0

	"""

		INITIALISATION METHODS

	"""

	def __init__(self, name, Mr, conc_mol_per_cell=None, conc_mol_per_l=None, thermo=True, gamma=1., T=298.15, P=101325.):
		self.name = name
		self.Mr = Mr

		if thermo: #pass Thermo as False to update thermochemical parameters yourself
			self.GetThermoParams(T, P)

		self.thermo = thermo
		self.conc_mol_per_cell = conc_mol_per_cell
		self.conc_mol_per_l = conc_mol_per_l
		self.gamma = gamma

	def GetThermoParams(self, T, P):

		db = Database(this_dir+"/../../../data/synthesis/supcrt07-organics_AABB.xml")
		thermodynamics = Thermo(db)
		self.std_formation_gibbs = thermodynamics.standardPartialMolarGibbsEnergy(T, P, self.name).val
		AABB = thermodynamics.standardPartialMolarGibbsEnergy(T, P, 'AABB').val
		# self.std_formation_gibbs = thermodynamics.standardPartialMolarGibbsEnergy(T=T, P=P, species=self.name).val
		# AABB = thermodynamics.standardPartialMolarGibbsEnergy(T=T, P=P, species='AABB').val
		self.std_formation_R = self.std_formation_gibbs - AABB
