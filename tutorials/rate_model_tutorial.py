import sys
# sys.path.append('..')
import NutMEG as nm
import NutMEG.core as nmc
import NutMEG.models as nmm


import numpy as np
import matplotlib.pyplot as plt


R = nmc.Reactor()
c = nmm.base_rate_models.ConstantRate(10.)
a = nmm.base_rate_models.Arrhenius(A=1000., Ea=13462)
q10 = nmm.base_rate_models.Q10(10., Q=2)


RMs = { 'Const' : (c, [], '-'),
  'Arrhenius' : (a, [], 'dotted'),
  'Q10' : (q10, [], 'dashed'),
}

Ts = np.linspace(273, 333, num=500)


for _T in Ts:
    R.T =  _T
    for k, v in RMs.items():
        v[1].append(v[0].compute(None, R))

fig, ax = plt.subplots(figsize=(6,4))
for k, v in RMs.items():
    ax.plot(Ts, v[1], label=k, ls=v[2], lw=3)
plt.xlim(270, 335)
ax.legend()
ax.set_ylabel('Rate [mol / s]')
ax.set_xlabel('Temperature [K]')


plt.savefig('tutorials_figures/base_rates.png')
plt.close()
