from .arrhenius import Arrhenius
from .base_rate_model import BaseRateModel
from .constant import Constant as ConstantRate
from .first_order_chemical import FirstOrderChemical

__all__ = ['BaseRateModel', 'Arrhenius', 'ConsantRate', 'FirstOrderChemical']


def find_brm(key):
    if key == 'BaseRateModel':
        return BaseRateModel
    if key == 'Arrhenius':
        return Arrhenius
    if key == 'ConstantRate':
        return ConstantRate
    if key == 'FirstOrderChemical':
        return FirstOrderChemical
