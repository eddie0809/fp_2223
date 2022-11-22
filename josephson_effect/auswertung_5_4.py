#!/usr/bin/env python3

import os
import glob
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import sinc
from scipy.constants import elementary_charge, k


def fit_sinc(x, a, b):
    return a * np.abs(sinc(b*x)) # sinc is defined with pi in scipy.special


def fit_tanh(x,a):
    return a * np.tanh(a*0.5*1.76/x)


def implicit_func(y,x):
    return np.tanh(y/x)


def lin_fit(x,a,b):
    return a*x+b


plt.rcParams.update({'text.usetex': True,})
plt.rcParams['text.latex.preamble'] = r'\usepackage[utf8]{inputenc} \usepackage[T1]{fontenc} \usepackage{amsmath} \usepackage{siunitx}'

cmap = mpl.cm.get_cmap('inferno') #['viridis', 'plasma', 'inferno', 'magma', 'cividis']
cmap_2 = mpl.cm.get_cmap('cividis')

cm = 1/2.54 # matplotlib is stupid and uses inches instead of cm.

files_5_3 = [file for file in glob.glob("files/5,3_*.txt")]
files_5_3.sort()
temps = [4.2, 4.82, 0.5*(5.68+5.96), 0.5*(6.47+6.55), 0.5*(7.01+7.51), 0.5*(7.55+7.99)]
i=0
i_c=[]
i_c_real=[]

#fig, ax = plt.subplots(figsize=(14*cm,14*0.75*cm))
for file in files_5_3[5:]:
    workinonit = (np.loadtxt(file, delimiter='\t', skiprows=1))*1e3
    newdiff = np.array([u - workinonit[2001,0] for u in workinonit[:,0]])
    sd_diff = 0.5*np.std(np.diff(newdiff))
    #ax.set_xlabel(r'$U$ [\si{\milli\volt}]')
    #ax.set_ylabel(r'$I_c$ [\si{\milli\ampere}]')
    #ax.set_xlim(1.8,2.2)
    #ax.set_ylim(.8,1.3)
    #ax.plot(newdiff, workinonit[:,1], label=r'$I_c(U)$, at $T$ = \SI{'+str(temps[i])+'}{\kelvin}', color=cmap((i+1)/8))
    #popt1, pcov1 = curve_fit(lin_fit, newdiff[:800], workinonit[:800,1])
    #popt2, pcov2 = curve_fit(lin_fit, newdiff[1000:1300], workinonit[1000:1300,1])
    #popt3, pcov3 = curve_fit(lin_fit, newdiff[830:900], workinonit[830:900,1])
    #ax.plot(newdiff, lin_fit(newdiff, *popt3))
    #ax.plot(newdiff, lin_fit(newdiff, *popt1))
    #ax.plot(newdiff, lin_fit(newdiff, *popt2))
    #print(popt1,popt2,popt3)
    #intersect = [
    #np.argwhere(np.diff(np.sign((-1*lin_fit(newdiff, *popt3) + lin_fit(newdiff, *popt1))))).flatten(),
    #np.argwhere(np.diff(np.sign((-1*lin_fit(newdiff, *popt3) + lin_fit(newdiff, *popt2))))).flatten()
    #]
    #print(0.5*(newdiff[intersect[0][0]]+newdiff[intersect[1][0]]))
#ax.legend(loc='best', ncol=3)
#fig.subplots_adjust()
#plt.savefig("plots/Ic_of_U_temp_dependence.pdf")
#plt.show()
#i_c_real[1] = 0.082 # both values read out manually from file
#i_c_real[2] = 0.076 # - " -

# für T = 4.20: 2.906     meV
# für T = 4.82: 2.8440938 meV
# für T = 5.68: 2.6779494 meV
# für T = 6.47: 2.4800684 meV
# für T = 7.01: 2.2298165 meV
# für T = 7.55: 1.9897503 meV
gaps = [2.906, 2.8440938, 2.6779494, 2.4800684, 2.2298165, 1.9897503]
y = np.linspace(0.01,1,1000, endpoint=True)
x = np.linspace(0.01,1,1000, endpoint=True)
new_y = []
y_sol=[]
for x in x:
    new_y.append(implicit_func(y,x))
#print(new_y)
#np.argwhere(np.diff(np.sign(f - g))).flatten()
for y_of_x in new_y[:-1]:
    y_sol.append(.001*(1+(np.argwhere(np.diff(np.sign((-1*np.array(y_of_x) + y))).flatten()))[0][0]))
#print(y_sol)
y_sol = np.array(y_sol)

fig, (ax1, ax) = plt.subplots(2,1,figsize=(14*cm, 63*cm/4))
i=0
ax1.plot(temps, gaps, '.', color = cmap(2/8), label=r'$2\Delta(T)$')
ax1.set_xlabel(r'$T$ [\si{\kelvin}]')
ax1.set_ylabel(r'$2\Delta(T)$ [\si{\milli\electronvolt}]')
ax.set_xlabel(r'$T/T_c$')
ax.set_ylabel(r'$\Delta(T)/\Delta_0$')
ax.plot(np.array(temps)*9.7**-1, np.array(gaps)*gaps[0]**-1, '.', color = cmap((5+1)/8), label=r'$\Delta(T)/\Delta_0$, with $T_c=\SI{9.7}{\kelvin}$')
ax.plot(np.array(temps)*9.2**-1, np.array(gaps)*gaps[0]**-1, '.', color = cmap((3+1)/8), label=r'$\Delta(T)/\Delta_0$, with $T_c=\SI{9.2}{\kelvin}$')
    #ax.plot(temp/9.2, i_c_real[i]/i_c_real[0], '.', color = cmap((i+1)/8))
    #ax.plot(temp/12, i_c_real[i]/i_c_real[0], '.', color = cmap((i+1)/8))
#ax.plot(T, fit_tanh(T, y), color='black')
ax.axhline(0,color='grey')
ax1.axhline(0,color='grey')
ax.set_xlim(.3, 1)
y=np.linspace(0.01,1,999)
ax.plot(y, y_sol, color=cmap(1), label = r'BCS-Theory', lw=1)
fig.tight_layout()
fig.subplots_adjust(bottom=.18)
fig.legend(loc=8, ncol=2)
plt.savefig("plots/energygap_of_T_done_ithink.pdf")
plt.show()