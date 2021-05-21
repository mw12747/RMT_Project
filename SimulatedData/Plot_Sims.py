import json
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt, pi


def plot_sims(hard, dim, rate, mean_degree, eta=2):
    if hard:
        fname = f'/Users/michaelwilsher/Documents/GitHub/RMT_Project/SimulatedData/HardRGG/' \
                f'Dim_{dim}/rate_{rate}/eta_inf.txt'
    else:
        fname = f'/Users/michaelwilsher/Documents/GitHub/RMT_Project/SimulatedData/SoftRGG/' \
                f'Dim_{dim}/rate_{rate}/eta_{eta}.txt'

    with open(fname, "r") as my_file:
        data = json.load(my_file)

    bins = np.arange(-40.05, 40.05, 0.1)

    if dim == 1:
        if hard == True:
            rad = mean_degree / (2 * rate)

        else:
            if eta == 1:
                rad = mean_degree / (2 * rate)
            else:
                rad = mean_degree / (rate * sqrt(pi))

    elif dim == 2:
        if hard == True:
            rad = sqrt(mean_degree / (pi * rate))

        else:
            if eta == 1:
                rad = sqrt(mean_degree / (2 * pi * rate))
            elif eta == 2:
                rad = sqrt(mean_degree / (pi * rate))

    key = str(rad)
    n_sims = len(data[key])
    print('nsims = ', n_sims)
    eig = []

    for i in range(n_sims):
        eig = eig + data[key][i]
    print(len(eig))

    plt.hist(eig, bins=len(eig), density=False, cumulative=True, label='CDF',
             histtype='step', alpha=0.8, color='k')

    # plt.hist(eig, bins=bins, density=True, label=f'mean degree = {mean_degree}')
    # if hard:
    #     plt.title(f'Hard RGG for d = {dim}, rate = {rate}')
    # else:
    #     plt.title(f'Soft RGG for d = {dim}, rate = {rate}, $\eta$ = {eta}')


    plt.legend()
    plt.show()


rate = 1e3
mean_degree = 10

plot_sims(hard=True, dim=1, rate=rate, mean_degree=mean_degree, eta=1)
# plot_sims(hard=False, dim=1, rate=rate, mean_degree=mean_degree, eta=2)
# plot_sims(hard=False, dim=1, rate=rate, mean_degree=mean_degree, eta=1)
