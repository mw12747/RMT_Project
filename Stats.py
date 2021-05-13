import pandas as pd
import matplotlib.pyplot as plt
from RMT_Project.RMTstats import Spacings, Spacings2
import numpy as np
import scipy as sp

hard = False
dim = 1
rate = 1e3
rad = 0.04
eta = 2

if hard:
    fname = '/Users/michaelwilsher/PycharmProjects/1DSoftRGG/' \
        'RMT_Project/SimulatedData/HardRGG/Dim_%s/rate_%s/rad_%s/eta_inf.csv' % (dim, rate, rad)
    fname_readme = '/Users/michaelwilsher/PycharmProjects/1DSoftRGG/RMT_Project/SimulatedData/HardRGG/' \
                   'Dim_%s/rate_%s/rad_%s/readme.txt' % (dim, rate, rad)
else:
    fname = '/Users/michaelwilsher/PycharmProjects/1DSoftRGG/' \
            'RMT_Project/SimulatedData/SoftRGG/Dim_%s/rate_%s/rad_%s/eta_%s.csv' % (dim, rate, rad, eta)
    fname_readme = '/Users/michaelwilsher/PycharmProjects/1DSoftRGG/RMT_Project/SimulatedData/SoftRGG/' \
                   'Dim_%s/rate_%s/rad_%s/readme_eta_%s.txt' % (dim, rate, rad, eta)

file = open(fname_readme, "r")
for line in file.readlines():
    if 'nsims' in line:
        nsims_temp = line[8:]
        nsims = float(nsims_temp)

df = pd.read_csv(fname, index_col=False)

data = {}
for key in df.keys():
    data[key] = df[key].tolist()

bins = np.arange(-40.05, 40.05, 0.1)

NNSD = np.array([]) # array to store the nearest spacings
nNNSD = np.array([]) # array to store next nearest spacings

for key in [str(rad)]:
    eig_temp = data[key]
    eig = [x for x in eig_temp if np.isnan(x) == False]
    sort_eig = np.sort(eig)
    unfold = sp.interpolate.interp1d(sort_eig, np.linspace(1.0, rate, sort_eig.size), bounds_error=False,
                                     fill_value=np.nan)
    NNSD = np.append(NNSD, Spacings(sort_eig))  # Calculate the NNSD
    nNNSD = np.append(nNNSD, Spacings2(sort_eig))  # Calculate the nNNSD
