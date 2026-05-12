from .bioenergetic_growth_model import BioenergeticGrowthModel
from .classic_growth_model import ClassicGrowthModel
from .growth_model import GrowthModel

__all__ = ['GrowthModel', 'ClassicGrowthModel', 'BioenergeticGrowthModel']

def find_gm(key):
    if key == 'GrowthModel':
        return GrowthModel
    if key == 'ClassicGrowthModel':
        return ClassicGrowthModel
    if key == 'BioenergeticGrowthModel':
        return BioenergeticGrowthModel
