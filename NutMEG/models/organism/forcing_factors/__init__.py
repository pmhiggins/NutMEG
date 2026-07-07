from .bioenergetic import Bioenergetic
from .forcing_factor import ForcingFactor
from .monod import Monod
from .combo_add import ComboAdd

__all__ = ['ForcingFactor', 'Bioenergetic', 'Monod']

def find_ff(key):
    if key == 'ForcingFactor':
        return ForcingFactor
    if key == 'Bioenergetic':
        return Bioenergetic
    if key == 'Monod':
        return Monod
    if key == 'ComboAdd':
        return ComboAdd
