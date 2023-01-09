#!/usr/bin/env python3

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from numpy import log as ln

plt.rcParams.update({'text.usetex': True,})
plt.rcParams['text.latex.preamble'] = r'\usepackage[utf8]{inputenc} \usepackage[T1]{fontenc} \usepackage{amsmath} \usepackage{siunitx}'

cmap = mpl.cm.get_cmap('plasma')


def paschens_law_fit(p, B, A, g):
    return B*p*20/(ln(A*p*20)-ln(ln(1+g))+ln(p*20))


p = (np.loadtxt("data/paschen_argon", skiprows=3)[:,1])*100
U = np.loadtxt("data/paschen_argon", skiprows=3)[:,0]
d = 0.2

pd = p*d*100

popt, pcov = curve_fit(paschens_law_fit, p, U)#, p0=[])

print(popt)
sd = np.square(np.diag(pcov))
print(sd)

p_lin = np.linspace(50, 300, 1000)

plt.plot(p, U, '.', color=cmap(1/3))
plt.plot(p_lin, paschens_law_fit(p_lin, *popt), color=cmap(2/3))
plt.xlabel(r'$p$ in pa')
plt.ylabel(r'$U$ in V')
plt.show()
