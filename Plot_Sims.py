import json
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt, pi

def mean_degree2rad(mean_degree, hard, rate, dim, eta=2):
    if dim == 1:
        if hard == True:
            rad = mean_degree / (2 * rate)
        else:
            if eta == 1:
                rad = mean_degree / (2 * rate)
            elif eta == 2:
                rad = mean_degree / (rate * sqrt(pi))
    elif dim == 2:
        if hard == True:
            rad = sqrt(mean_degree / (pi * rate))
        else:
            if eta == 1:
                rad = sqrt(mean_degree / (2 * pi * rate))
            elif eta == 2:
                rad = sqrt(mean_degree / (pi * rate))
    return rad

def rad2mean_degree(rad, hard, rate, dim, eta=2):
    if dim == 1:
        if hard:
            mean_degree = 2 * rate * rad
        else:
            if eta == 1:
                mean_degree = 2 * rate * rad
            elif eta == 2:
                mean_degree = rate * sqrt(pi) * rad
    elif dim == 2:
        if hard:
            mean_degree = pi * rate * rad ** 2
        else:
            if eta == 1:
                mean_degree = 2 * pi * rate * rad ** 2
            elif eta == 2:
                mean_degree = pi * rate * rad ** 2
    return mean_degree

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

    rad = mean_degree2rad(mean_degree, hard, rate, dim, eta)

    # rad_list_temp = sorted(data.keys())
    # rad_list = [float(x) for x in rad_list_temp]
    # mean_degrees_list = [rad2mean_degree(x, hard, rate, dim, eta) for x in rad_list]
    # print(mean_degrees_list)

    key = str(rad)
    n_sims = len(data[key])
    print('nsims = ', n_sims)
    eig = []

    for i in range(n_sims):
        eig = eig + data[key][i]
    # print(len(eig))

    plt.hist(eig, bins=bins, density=True, label='rad = %s' % key)
    if hard:
        plt.title(f'Hard RGG for d = {dim}, rate = {rate}')
    else:
        plt.title(f'Soft RGG for d = {dim}, rate = {rate}, $\eta$ = {eta}')


    plt.legend()
    plt.show()


# plot_sims(True, 2, 1e3, 125, 1)
# plot_sims(False, 2, 1e3, 125, 2)
plot_sims(False, 2, 1e3, 125, 1)
#
print(rad2mean_degree(0.1, True, 1e3, 2))