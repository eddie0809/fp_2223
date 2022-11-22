#!/usr/bin/env python3

import os
import glob
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import sinc


def fit_sinc(x, a, b):
    return a * np.abs(sinc(b*x/np.pi)) # sinc is defined with pi in scipy.special

plt.rcParams.update({'text.usetex': True,})
plt.rcParams['text.latex.preamble'] = r'\usepackage[utf8]{inputenc} \usepackage[T1]{fontenc} \usepackage{amsmath} \usepackage{siunitx}'

cmap = mpl.cm.get_cmap('inferno') #['viridis', 'plasma', 'inferno', 'magma', 'cividis']
cmap_2 = mpl.cm.get_cmap('cividis')

cm = 1/2.54 # matplotlib is stupid and uses inches instead of cm.

path = r'~/Documents/fp_2223/josephson_effect/files/'
data = []

i_v_b_table = np.loadtxt("files/I-V-B-table.csv", delimiter='\t', skiprows=1)

B_values = np.append(i_v_b_table[:,2], (i_v_b_table[1:,2])*-1)
#B_values = np.delete(B_values, 1)
print(B_values)

folder_name = 'files'
file_type = 'txt'

i_c = []
i_c_real = []
data_5_2 = {}
stepsize = 1e-6

files_5_2 = [file for file in glob.glob(folder_name + "/5,2_*."+file_type)]
files_5_2.sort()

for file in files_5_2:
    workinonit = np.loadtxt(file, delimiter='\t', skiprows=1)*1e3
    newdiff = [u - workinonit[201,0] for u in workinonit[:,0]]
    sd_diff = .5*np.std(np.diff(newdiff))**2
    print(sd_diff)
    print()
    i_c.append(np.argwhere(np.array(newdiff)<sd_diff)[0][0])
    i_c.append(np.argwhere(np.flip(newdiff)<sd_diff)[-1][0])
    i_c_real.append(np.abs(workinonit[i_c[-2],1]-workinonit[i_c[-1],1])/2)
    data_5_2[str(file)] = workinonit

err = 1e-5*np.ones(len(i_c_real))
x = np.linspace(B_values.min(),B_values.max(),300)

popt, pcov = curve_fit(fit_sinc, B_values, i_c_real, p0=(0.086, 2.6)) #[8.29622403e-05 8.76209175e-01]

print("Amplitude a = " + str(popt[0])+"\n")
print("width b = "+str(popt[1])+"\n")
print("standard deviations [a,b]: "+ str(np.square(np.diag(pcov)))+'\n')
print("critical currents: " + str(i_c_real) + '\n')
#table = pd.DataFrame(i_v_b_table)
#print(table)
#print(minima, diff_u)
#print()
#print(i_c)
#data_5_2 = np.array([pd.read_table(f, delimiter='\t', skiprows=1, names=str(f)) for f in glob.glob(folder_name + "/5,2_*."+file_type)])
#print(data_5_2)

if __name__=='__main__':
    fig, ax = plt.subplots(1,1, figsize=(14*cm,9*7*cm/8))

    ax.set_xlabel(r'$B$ in $\si{\milli\tesla}$')
    ax.set_ylabel(r'$I_c$ in $\si{\milli\ampere}$')
    #
    ax.plot(x, fit_sinc(x, *popt), color = cmap(1/8), label = r'$f(B)$', alpha=0.8)
    ax.plot(B_values, i_c_real, '+', color = cmap(9/16), label = r'$I_c(B)$', alpha=1)
    #plt.show()
    ##ax.plot(data_5_2['files/5,2_M01_0A.txt'][:,0], data_5_2['files/5,2_M01_0A.txt'][:,1], '.', label = r'$I_c$', color=cmap(1/2))
    #ax2.plot(data_5_2[str('files/5,2_M23_-1,0A.txt')][:,0], data_5_2[str('files/5,2_M23_-1,0A.txt')][:,1], '.', markersize = 1, color = cmap(1/2))
    #ax2.axhline(data_5_2[str('files/5,2_M23_-1,0A.txt')][i_c[-5],1], color = 'grey')
    #ax2.axhline(data_5_2[str('files/5,2_M23_-1,0A.txt')][i_c[-6],1], color = 'grey')
    ax.legend(loc='best')
    #plt.savefig('test_Ic_of_B.pdf')
    fig.subplots_adjust()
    fig.tight_layout()
    plt.show()