#!/usr/bin/env python3

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from numpy import log as ln
from numpy import exp
from scipy.constants import epsilon_0, e, k, m_e
from lmfit import Model

plt.rcParams.update({'text.usetex': True,})
plt.rcParams['text.latex.preamble'] = r'\usepackage[utf8]{inputenc} \usepackage[T1]{fontenc} \usepackage{amsmath} \usepackage{siunitx}'

cmap = mpl.cm.get_cmap('plasma')


def lin_fit(U, a, b):
    return a*U+b


def n_e(I_sat_e, T_e, S = 1.3e-6, m=m_e):
    return I_sat_e*np.sqrt(2*np.pi*m)/(e*S*np.sqrt(T_e))


U = np.loadtxt("data/Argon2mbar0R.txt")[:,1]
I = ln(np.loadtxt("data/Argon2mbar0R.txt")[:,2]*-1000+0.011)

I_sorted = np.sort(I)
U_sorted = np.sort(U)

sign = np.sign(I_sorted)
fl_ind = np.nonzero(np.diff(sign))[0][0]
fit_ind = np.argmax(np.diff(U_sorted))
U_float = U_sorted[fl_ind]
print(U_sorted[fit_ind-10:fit_ind])
popt, pcov = curve_fit(lin_fit, U_sorted[fit_ind-32:fit_ind-17], I_sorted[fit_ind-32:fit_ind-17])#, p0=[1/6, 23, e*U_float])

print(popt)
print(np.square(np.diag(pcov)))
T_e = popt[0]**-1
n = n_e(exp(popt[1]-e*U_float*popt[0])*1e3, T_e)

print(f"T_e = {T_e}, n_e = {n}")


U_lin  = np.linspace(U_sorted[0], U_sorted[fit_ind], 200)
fig = plt.figure(figsize=(14.4/2.54, 63/(8*2.54)))
plt.plot(U, I, '.', color=cmap(1/3))
plt.plot(U_lin, lin_fit(U_lin, *popt), color=cmap(2/3))
#plt.yscale('log')
plt.ylabel(r'$-(I-I_\mathrm{i,sat})$ in mA')
plt.xlabel(r'$U$ in V')
fig.tight_layout()
#plt.savefig("plots/Argon_charakteristik_example.pdf")
plt.show()