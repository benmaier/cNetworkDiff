import numpy as np
from math import log10
import pylab as pl
import seaborn as sns
import cNetworkDiffusion as diff
import progressbar

sns.set_style("whitegrid")

B = 8
L = 3
k = 7
N_xi = 50 
xis = np.logspace(log10(0.25),log10(B),N_xi)
coverage_ratio = 1.0

N_meas = 10

seed = 23556

mmfpt_pow = np.zeros(N_xi)
mcovt_pow = np.zeros(N_xi)

mmfpt_mhr = np.zeros(N_xi)
mcovt_mhr = np.zeros(N_xi)

bar = progressbar.ProgressBar(max_value=N_xi-1)
for ix,xi in enumerate(xis):

    for meas in range(N_meas):

        m,c = diff.mmfpt_and_mean_coverage_time_meanfield_MHRN(B,L,k,xi,coverage_ratio,seed)

        mmfpt_mhr[ix] += m
        mcovt_mhr[ix] += c

        seed += 1


        alpha = np.log(xi)/np.log(B)-1 

        #m,c = diff.mmfpt_and_mean_coverage_time_power_law_k(B**L,alpha,k,coverage_ratio,seed)
        m,c = diff.mmfpt_and_mean_coverage_time_power_law(B**L,alpha,coverage_ratio,seed)

        mmfpt_pow[ix] += m
        mcovt_pow[ix] += c

        seed += 1
    bar.update(ix)


mmfpt_mhr /= N_meas
mcovt_mhr /= N_meas
mmfpt_pow /= N_meas
mcovt_pow /= N_meas

fig, ax = pl.subplots(1,2,figsize=(12,6))

#pl.yscale("log")
ax[0].plot(np.log(xis)/np.log(B),mmfpt_mhr,'k-',label='MHRN pmf')
ax[0].plot(np.log(xis)/np.log(B),mmfpt_pow,'k--',label=r'$|\Delta x|^{\mu-1}$ pmf with $\mu=\log\xi/\log B$')
ax[0].set_ylabel(r"mean global mean first passage time [steps]")

ax[1].plot(np.log(xis)/np.log(B),mcovt_mhr,'k-',label='MHRN pmf')
ax[1].plot(np.log(xis)/np.log(B),mcovt_pow,'k--',label=r'$|\Delta x|^{\mu-1}$ pmf with $\mu=\log\xi/\log B$')
ax[1].set_ylabel(r"mean coverage time [steps]")

for iax in range(2):
    ax[iax].set_xlim(np.log(xis)[[0,-1]]/np.log(B))
    ax[iax].legend(loc='best')
    ax[iax].set_xlabel(r"hierarchical structure parameter $\log\xi/\log B$")

#pl.plot(np.log(xis)/np.log(10),mcovt)
fig.savefig("mgmfpt_simulation_MHRN_pow.pdf")
pl.show()

