#!/usr/bin/env python3

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from numpy import log as ln
from numpy import exp
from scipy.constants import epsilon_0, e, k, m_e


plt.rcParams.update({'text.usetex': True,})
plt.rcParams['text.latex.preamble'] = r'\usepackage[utf8]{inputenc} \usepackage[T1]{fontenc} \usepackage{amsmath} \usepackage{siunitx}'

cmap = mpl.cm.get_cmap('plasma')

def lin_fit(U, a, b):
    return a*U+b


def I_sat(S, m, n, T):
    return n * e * S * np.sqrt(T/m)


def characteristics_fit(U, T_e, n, phi_pl): # fit between -5 and 20 volts
    S=1.3e-6        # convert mm^2 to m^2
    # T_e is in eV in the formula, we need to convert to Kelvin at the end
    return -(2*np.pi)**(-1/2)*I_sat(S, m_e, n, T_e) * exp(-e*(phi_pl-U)/T_e)

U = np.loadtxt("data/ArgonCharakteristik0,4mbar(20mA500V).txt")[:,1]
I = np.loadtxt("data/ArgonCharakteristik0,4mbar(20mA500V).txt")[:,2]*-1000
wheretofit = np.argwhere(U>-5)
#popt, pcov = curve_fit(characteristics_fit, np.asarray(U>-5 and U<20).nonzero(), I)
#print(np.asarray(U>-5 and U<20).nonzero())
print(U_fit)
plt.plot(U, I, '.', color=cmap(1/3))
#plt.plot(p_lin, paschens_law_fit(p_lin, *popt), color=cmap(2/3))
plt.ylabel(r'$I$ in $\si{\milli\ampere}$')
plt.xlabel(r'$U$ in V')
plt.show()
