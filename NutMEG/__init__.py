# from .reactor import reactor
# from .environment import environment
# from .colony import colony
# from .oranism import organism
# from .reaction import reaction
from .ecosystem import ecosystem
from .ecosystem_dbhelper import db_helper
from .environment import environment



#culture stuff
from .culture.base_organism.base_organism import base_organism
from .culture.base_organism.base_organism_dbhelper import bodb_helper
from .culture.organism.organism import organism
from .culture.horde.horde import horde
from .culture.colony.colony import colony
from .culture.culture import culture

from .reactor.reactor import reactor

from .culture.base_organism.synthesis.BioMolecule import BioMolecule


# from .util.loggersetup import loggersetup as logset
# loggers = logset(__name__, level='DEBUG')

# logger = logset.get_logger(__name__, level='DEBUG')
# import reaction as rxn
# #from reaction.reaction import reaction
# from reaction.reagent import reagent
# from environment.environment import environment
# from organism.organism import organism
# from organism.respirator import respirator
# from organism.CHNOPSexchanger import CHNOPSexchanger
# from reactor.reactor import reactor
# from reactor.reactormesh import reactormesh
# from colony.colony import colony
# from colony.colony_lite2 import colony_lite2
# from ecosystem import ecosystem
# from ecosystem_dbhelper import db_helper
#
# sys.path.append(this_dir+'/organism')
# sys.path.append(this_dir+'/reactor')
#from organism.saved_organisms import Methanogen2
#from reactor.saved_systems import Enceladus
