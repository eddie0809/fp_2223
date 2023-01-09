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


def I_sat(S, m, n, T):
    return n * e * S * np.sqrt(T/m)


def n_e(I_sat_e, T_e, S = 1.3e-6, m=m_e):
    return I_sat_e*np.sqrt(2*np.pi*m)/(e*S*np.sqrt(T_e))


def alpha(n, T_G, p_G, kb=k):
    return n*T_G*kb/p_G


def debye_length(T, n):
    return np.sqrt(epsilon_0 *k * T/(n*e**2))


def plasma_freq(n):
    return np.sqrt(n*e**2/(epsilon_0*m_e))


def characteristics_fit(x, T_e, a, phi_pl): # fit between -5 and 20 volts
    #S=1.3e-6         convert mm^2 to m^2
    # T_e is in eV in the formula, we need to convert to Kelvin at the end * I_sat(1.3e-6, 1, n, T_e)
    return a * exp(-e*(phi_pl-x)*T_e)


def char_fit(x, T, n, U_fl):
    return 0.61*n*e*1.3e-6*np.sqrt(T/9.458e-24) * (1-exp(e*(x-U_fl)/T))

#0.61*n*e*1.3e-6*np.sqrt(T/9.458e-24)

U = np.loadtxt("data/HeliumCharakteristik7mbar.txt")[:,1]
I = ln(np.loadtxt("data/HeliumCharakteristik7mbar.txt")[:,2]*-1000+0.048)

I_sorted = np.sort(I)
U_sorted = np.sort(U)

I_sat_i = I_sorted[0]
sign = np.sign(I_sorted)
fl_ind = np.nonzero(np.diff(sign))[0][0]
fit_ind = np.argmax(np.diff(U_sorted))
U_float = U_sorted[fl_ind]
print(U_sorted[fit_ind-10:fit_ind])
popt, pcov = curve_fit(lin_fit, U_sorted[fit_ind-17:fit_ind-1], I_sorted[fit_ind-17:fit_ind-1])#, p0=[1/6, 23, e*U_float])

print(popt)
print(np.square(np.diag(pcov)))
T_e = popt[0]**-1
n = n_e(exp(popt[1]-e*U_float*popt[0])*1e3, T_e)
alph = alpha(n, 293, 700)
lambda_d = debye_length(T_e, n)
omega_p = plasma_freq(n)
print(f"T_e = {T_e}, n_e = {n}, alpha = {alph}, lambda_d = {lambda_d}, omega_p = {omega_p}")
#gmodel = Model(char_fit)
#result = gmodel.fit(I_sorted[U_sorted<0], x=U_sorted[U_sorted<0], nan_policy='propagate')#, T=8, n, U_fl = U_float)
#print(result.fit_report())

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
