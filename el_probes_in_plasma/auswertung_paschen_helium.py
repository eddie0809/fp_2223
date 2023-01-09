#!/usr/bin/env python3

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from numpy import log as ln

plt.rcParams.update({'text.usetex': True,})
plt.rcParams['text.latex.preamble'] = r'\usepackage[utf8]{inputenc} \usepackage[T1]{fontenc} \usepackage{amsmath} \usepackage{siunitx}'

cmap = mpl.cm.get_cmap('plasma')


def paschens_law_fit(pd, B, A, g):
    return B*pd/(ln(A*pd)-ln(ln(1+g))+ln(pd))


def paschens_law_fit_less_fitparam(pd, B, A):
    return B*pd/(ln(A*pd))


def paschens_law_fit_no_ln(pd, B, A, g):
    return B*pd/(ln(A*pd)-ln(ln(1+g)))


p = (np.loadtxt("data/paschen_helium", skiprows=1)[:,1])*100
U = np.loadtxt("data/paschen_helium", skiprows=1)[:,0]
d = 0.2

pd = p*d

# stainless steel electrodes, (for A, B, g)

#popt, pcov = curve_fit(paschens_law_fit, pd[0:], U[0:])#, p0=[])
popt1, pcov1 = curve_fit(paschens_law_fit_no_ln, pd[1:], U[1:])

print(popt1)
sd1 = np.square(np.diag(pcov1))
print(sd1)

p_lin = np.linspace(p[1]*d, p[-1]*d, 1000)

plt.plot(pd, U, '.', color=cmap(1/4))
#plt.plot(p_lin, paschens_law_fit(p_lin, *popt), color=cmap(1/2))
plt.plot(p_lin, paschens_law_fit_no_ln(p_lin, *popt1), color=cmap(3/4))
plt.xlabel(r'$p$ in pa')
plt.ylabel(r'$U$ in V')
plt.show()
