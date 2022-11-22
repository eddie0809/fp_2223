#!/usr/bin/env python3

import glob
import pandas as pd
import numpy as np

data = np.loadtxt("files/I-V-B-table.csv", skiprows=1)
data_pd = pd.read_csv("files/I-V-B-table.csv", delimiter='\t')#, header=1)

#print(data, data_pd)
#heft_werte = pd.DataFrame(data_pd)
print(data_pd)
print()
table = data_pd.to_latex(index=False, header=True, column_format='ccc')
print(table)