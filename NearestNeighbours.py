import json
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt, pi
import scipy as sp
import scipy.interpolate


def data_import(hard, dim, rate, mean_degree, eta=2):
    """
    Imports data from .txt file as a dictionary
    :param hard:
    :param dim:
    :param rate:
    :param mean_degree:
    :param eta:
    :return: Dictionary containing simulation results
    """
    if hard:
        fname = f'/Users/michaelwilsher/Documents/GitHub/RMT_Project/SimulatedData/HardRGG/' \
                f'Dim_{dim}/rate_{rate}/eta_inf.txt'
    else:
        fname = f'/Users/michaelwilsher/Documents/GitHub/RMT_Project/SimulatedData/SoftRGG/' \
                f'Dim_{dim}/rate_{rate}/eta_{eta}.txt'

    with open(fname, "r") as my_file:
        data = json.load(my_file)

    return data


def mean_degree2rad(mean_degree, hard, rate, eta=2):
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

#
def NNSD(data, mean_degree, hard, rate, eta=2):
    rad = str(mean_degree2rad(mean_degree, hard, rate, eta))
    list_len = len(data[rad])
    dat_temp = [sorted(data[rad][i]) for i in range(list_len)]
    NN = []
    unfold_list = []
    sum = 0

    for i in range(list_len):
        len_sub_list = len(dat_temp[i])
        sum += len_sub_list
        unfold = sp.interpolate.interp1d(dat_temp[i], np.linspace(1.0, len_sub_list, len_sub_list), bounds_error=False,
                                         fill_value=np.nan)
        dat = sorted(unfold(dat_temp[i]))

        y = 1000.0 * (np.arange(1, len_sub_list + 1) / len_sub_list)

        plt.step(dat_temp[i], y, color='red')

        plt.show()

        unfold_list.append(dat)
        if i == 40:
            print("NNSD_temp = ",dat)
        # print(dat)
        for x in range(len_sub_list - 1):
            NN.append(dat[x + 1] - dat[x])
    # print(NN)
    return NN, unfold_list

# def NNSD(data, mean_degree, hard, rate, eta=2):
#     rad = str(mean_degree2rad(mean_degree, hard, rate, eta))
#     list_len = len(data[rad])
#     dat = [sorted(data[rad][i]) for i in range(list_len)]
#     print("NNSD = ",dat[40])
#     NN = [dat[i][x + 1] - dat[i][x] for i in range(list_len) for x in range(len(dat[i]) - 1)]
#     # print(sorted(NN)[0:1000])
#     return NN

def nNNSD(data, mean_degree, hard, rate, eta=2):
    rad = str(mean_degree2rad(mean_degree, hard, rate, eta))
    list_len = len(data[rad])
    dat = [sorted(data[rad][i]) for i in range(list_len)]

    nNN = [(dat[i][x + 2] - dat[i][x]) / 2 for i in range(list_len) for x in range(len(dat[i]) - 2)]
    # print(nNN)
    return nNN

# print(mean_degree2rad(10, True, 1e3))

# sort_aspec = [1,2,3,4,5,5]
# unfold = sp.interpolate.interp1d(sort_aspec, np.linspace(1.0, 1e3, len(sort_aspec)),bounds_error=False, fill_value=np.nan)
# print(unfold)

hard = True
dim = 2
rate = 1e3
mean_degree = 125

# data = data_import(hard=hard, dim=dim, rate=rate, mean_degree=mean_degree)
#
# NN, unfold_list = NNSD(data, mean_degree, hard, rate)

data = data_import(hard=False, dim=dim, rate=rate, mean_degree=mean_degree, eta=1)

NN, unfold_list = NNSD(data, mean_degree, hard, rate, eta=1)
# NN_temp = NNSD_temp(data, mean_degree, hard, rate)

# print(NN[3])
# nNN = nNNSD(data, mean_degree, hard, rate)
#
# # print('NN = ', NN)
# # print('nNN = ', nNN)
#
# bins = list(np.arange(0,10,0.001))
#
# plt.figure()
# plt.hist(NN, bins=bins)
#
# plt.figure()
# plt.hist(nNN, bins=bins)
#
# plt.show()


