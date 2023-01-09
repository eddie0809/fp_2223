#!/usr/bin/env python3

import os

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from numpy import log as ln

plt.rcParams.update({'text.usetex': True,})
plt.rcParams['text.latex.preamble'] = r'\usepackage[utf8]{inputenc} \usepackage[T1]{fontenc} \usepackage{amsmath} \usepackage{siunitx}'

cmap = mpl.cm.get_cmap('plasma')
cm = 2.54**-1

def paschens_law_fit(p, B, A, g):
    return B*p*20/(ln(A*p*20)-ln(ln(1+g))+ln(p*20))



path = r'./data/'
files = []
for f in os.listdir(path):
    files.append(f)
files.sort()

p_ar = (np.loadtxt(path + files[-2], skiprows=3)[:,1])*100
U_ar = np.loadtxt(path + files[-2], skiprows=3)[:,0]
p_he = (np.loadtxt(path + files[-1], skiprows=1)[:,1])*100
U_he = np.loadtxt(path + files[-1], skiprows=1)[:,0]

d = 0.2
pd_ar = p_ar*d*100
pd_he = p_he*d*100

popt_ar, pcov_ar = curve_fit(paschens_law_fit, p_ar, U_ar)#, p0=[])
popt_he, pcov_he = curve_fit(paschens_law_fit, p_he[1:], U_he[1:])#, p0=[])

print(f"argon, [B, A, gamma] = {popt_ar},\nhelium, [B,A,gamma] = {popt_he}\n")

sd_ar = np.square(np.diag(pcov_ar))
sd_he = np.square(np.diag(pcov_he))
print(f"sd, argon = {sd_ar}, \nsd, helium = {sd_he}\n")

p_lin_ar = np.linspace(.50, 3.00, 1000)
p_lin_he = np.linspace(2.750, 10.00, 1000)

fig, (ax_ar, ax_he) = plt.subplots(2,1, figsize=(14.4*cm,63*cm*0.2))

ax_ar.plot(p_ar*0.01, U_ar, '.', color=cmap(1/3), label=r"Data, Argon")
ax_ar.plot(p_lin_ar, paschens_law_fit(p_lin_ar*100, *popt_ar), color=cmap(2/3), label = r"Paschen's law fit, Argon")
ax_he.plot(p_he*0.01, U_he, '.', color = cmap(1/3), label=r"Data, Helium")
ax_he.plot(p_lin_he, paschens_law_fit(p_lin_he*100, *popt_he), color = cmap(2/3), label = r"Paschen's law fit, Helium")
#ax_ar.set_xlabel(r'$p$ in pa')
ax_he.set_xlabel(r'$p$ in mbar')
ax_ar.set_ylabel(r'$U$ in V')
ax_he.set_ylabel(r'$U$ in V')
ax_ar.legend(loc='best')
ax_he.legend(loc='best')
fig.tight_layout()
plt.savefig("plots/paschens_law.pdf")
plt.show()
