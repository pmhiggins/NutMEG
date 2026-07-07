from .bioenergetic import Bioenergetic
from .forcing_factor import ForcingFactor
from .monod import Monod
from .upper_limit import UpperLimit
from .lower_limit import LowerLimit
from .combo_add import ComboAdd


__all__ = ['ForcingFactor', 'Bioenergetic', 'Monod', 'UpperLimit', 'LowerLimit', 'ComboAdd']

def find_ff(key):
    if key == 'ForcingFactor':
        return ForcingFactor
    if key == 'Bioenergetic':
        return Bioenergetic
    if key == 'Monod':
        return Monod
    if key == 'UpperLimit':
        return UpperLimit
    if key == 'LowerLimit':
        return LowerLimit
    if key == 'ComboAdd':
        return ComboAdd
