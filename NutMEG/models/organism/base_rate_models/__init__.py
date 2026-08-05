from .arrhenius import Arrhenius
from .base_rate_model import BaseRateModel
from .constant import Constant as ConstantRate
from .first_order_chemical import FirstOrderChemical
from .q10 import Q10


__all__ = ['BaseRateModel', 'Arrhenius', 'ConsantRate', 'FirstOrderChemical', 'Q10']


def find_brm(key):
    if key == 'BaseRateModel':
        return BaseRateModel
    if key == 'Arrhenius':
        return Arrhenius
    if key == 'ConstantRate':
        return ConstantRate
    if key == 'FirstOrderChemical':
        return FirstOrderChemical
    if key == 'Q10':
        return Q10
