import sys
# sys.path.append('..')
import NutMEG as nm
import NutMEG.core as nmc
import NutMEG.models as nmm


import numpy as np
import matplotlib.pyplot as plt

R = nmc.Reactor()
nmc.Reagent('H2S', R, amount=(0.001, 'molal'), thermo=False)
BP = nmm.forcing_factors.BiologicalPerformance('H2S', 0.075, 0.01, 0.081, 2.)
M = nmm.forcing_factors.Monod('H2S', 0.05)
LL = nmm.forcing_factors.LowerLimit('H2S', 0.01)
UL = nmm.forcing_factors.UpperLimit('H2S', 0.081)

FFs = { 'BiologicalPerformance' : (BP, [], '-'),
  'Monod' : (M, [], '-'),
  'LowerLimit' : (LL, [], 'dashed'),
  'UpperLimit' : (UL, [], 'dotted'),
}

Ts = np.linspace(0.0, 0.1, num=500)
Fs = []

for _T in Ts:
    R.composition['H2S'].molality =  _T
    for k, v in FFs.items():
        v[1].append(v[0].compute(None, R))

fig, ax = plt.subplots(figsize=(6,4))
for k, v in FFs.items():
    ax.plot(Ts, v[1], label=k, ls=v[2], lw=3)
# plt.xlim(270, 350)
ax.set_ylim(0, 1.2)
ax.legend()
ax.set_ylabel('Forcing factor')
ax.set_xlabel('[H2S]')


plt.savefig('tutorials_figures/rate_forcing.png')
plt.close()
