import sys
sys.path.append('..')
import NutMEG as nm
import NutMEG.core as nmc
import NutMEG.models as nmm


import numpy as np
import matplotlib.pyplot as plt


R = nmc.Reactor()
nmc.Reagent('N', R, amount=(0.001, 'molal'), thermo=False)

# assign a constant max rate.
c = nmm.base_rate_models.ConstantRate(10.)

# Make two forcing factors, one for T, one for N concentration:
BP = nmm.forcing_factors.BiologicalPerformance('T', 300, 275, 305, 2.)
M = nmm.forcing_factors.Monod('N', 0.05)

Agg = nmm.aggregators.Multiplicative()
Agg1 = nmm.aggregators.LeibigMinimum()

# initialise a BaseOrganism with these forcing factors
org = nmc.BaseOrganism(
  'organism',
  nmc.org.Metaboliser(None, base_rate=c, forcing_factors={'BioPerf':BP, 'Monod':M}))


fig, axs = plt.subplots(ncols=2, figsize=(10,4), sharey=True)

Ns = [0.001, 0.05, 0.1, 0.5] # N concentrations
Ts = np.linspace(273, 333, num=500) # temperature range

for i, A in enumerate([Agg, Agg1]):
    # Assign the organism's metabolic aggregator
    org.metabolism.aggregator = A
    for N in Ns:
        R.composition['N'].update_amount(R, mol=N) # set the N concentration
        rates = []
        for T in Ts:
            R.T = T # set the locale temperature
            # calculate the rate. this applies the aggregator.
            rates.append(org.metabolism.compute_rate(org, R))
        axs[i].plot(Ts, rates, label=f'N = {N} M')


for ax in axs:
    ax.legend()
    ax.set_ylabel('Rate [mol / s]')
    ax.set_xlabel('Temperature [K]')
axs[0].set_title('Multiplicative aggregator')
axs[1].set_title('Minimum aggregator')


plt.savefig('examples_figures/aggregators.png')
plt.close()
