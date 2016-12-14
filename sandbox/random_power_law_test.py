import numpy as np
from numpy.random import random
import pylab as pl
import seaborn as sns

N = 10
x0 = 1
x1 = N+1
x = np.arange(1,N+1)
bins = np.arange(1,N+2) - 0.5


alpha = -0.5
C = (alpha+1.)/(x1**(alpha+1)-x0**(alpha+1))

pmf = lambda x: x**alpha
C2 = 1./np.sum(pmf(x)) 

u = random(100000)

X = np.floor(( (alpha+1)/C2*u+x0**(alpha+1))**(1./(alpha+1)) )

pl.hist(X,bins=bins,normed=True)
pl.plot(x,C2*pmf(x))

pl.show()
