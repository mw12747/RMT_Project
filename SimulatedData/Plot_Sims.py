import json
import matplotlib.pyplot as plt
import numpy as np


def plot_sims(hard, dim, rate, rad, eta=2):
    if hard:
        fname = f'/Users/michaelwilsher/Documents/GitHub/RMT_Project/SimulatedData/HardRGG/' \
                f'Dim_{dim}/rate_{rate}/eta_inf.txt'
    else:
        fname = f'/Users/michaelwilsher/Documents/GitHub/RMT_Project/SimulatedData/SoftRGG/' \
                f'Dim_{dim}/rate_{rate}/eta_{eta}.txt'

    with open(fname, "r") as my_file:
        data = json.load(my_file)

    bins = np.arange(-40.05, 40.05, 0.1)

    key = str(rad)
    n_sims = len(data[key])
    print('nsims = ', n_sims)
    eig = []

    for i in range(n_sims):
        eig = eig + data[key][i]
    # print(len(eig))
    #
    # plt.hist(eig, bins=bins, density=True, label='rad = %s' % key)
    # if hard:
    #     plt.title(f'Hard RGG for d = {dim}, rate = {rate}')
    # else:
    #     plt.title(f'Soft RGG for d = {dim}, rate = {rate}, $\eta$ = {eta}')
    #
    #
    # plt.legend()
    # plt.show()


plot_sims(True, 1, 30, 0.016666666666666666, 1)
plot_sims(False, 1, 30, 0.01880631945159188, 2)
plot_sims(False, 1, 30, 0.016666666666666666, 1)
