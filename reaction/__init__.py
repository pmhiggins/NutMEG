"""
from reaction import reaction
from reaction import reagent
	# something similarly funky is needed for the special submodule.
from special import Tdep
from special.solutions import neutralsol
from special.solutions import electrolyte
from special.solutions import redox_half
from special.solutions import redox
"""



# for python 3
#"""
	# when reaction is imported, all we need to initialise the objects is reaction.reaction or reaction.reagent.
from .reaction import reaction
from .reagent import reagent
# import special

# #from reaction.reactive_system import reactive_system
# #from reaction.saved_systems.Enceladus import Enceladus
# 	# something similarly funky is needed for the special submodule.
from .special.Tdep import Tdep
from .special.solutions import neutralsol
from .special.solutions import electrolyte
from .special.solutions import redox_half
from .special.solutions import redox
#"""
