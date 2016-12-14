import numpy as np
from math import log10
import pylab as pl
import seaborn as sns
import cNetworkDiffusion as diff
import progressbar

B = 10
N_xi = 50
xis = np.logspace(log10(0.25),log10(10),N_xi)
N = 500
coverage_ratio = 1.0
k = 6;

seed = 2353456

mmfpt = np.zeros(N_xi)
mcovt = np.zeros(N_xi)

bar = progressbar.ProgressBar(max_value=N_xi-1)
for ix,xi in enumerate(xis):
    alpha = np.log(xi)/np.log(10) - 1  
    m,c = diff.mmfpt_and_mean_coverage_time_power_law_k(N,alpha,k,coverage_ratio,seed)
    mmfpt[ix] = m
    mcovt[ix] = c
    seed += 1
    bar.update(ix)

pl.plot(np.log(xis)/np.log(10),mmfpt)
#pl.plot(np.log(xis)/np.log(10),mcovt)
pl.show()

