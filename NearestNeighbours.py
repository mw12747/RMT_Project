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
    if hard == True:
        rad = mean_degree / (2 * rate)
    else:
        if eta == 1:
            rad = mean_degree / (2 * rate)
        else:
            rad = mean_degree / (rate * sqrt(pi))
    return rad

#
def NNSD_temp(data, mean_degree, hard, rate, eta=2):
    rad = str(mean_degree2rad(mean_degree, hard, rate, eta))
    list_len = len(data[rad])
    dat_temp = [sorted(data[rad][i]) for i in range(list_len)]
    NN = []
    sum = 0
    for i in range(list_len):
        len_sub_list = len(dat_temp[i])
        sum += len_sub_list
        unfold = sp.interpolate.interp1d(dat_temp[i], np.linspace(1.0, rate, len_sub_list), bounds_error=False,
                                         fill_value=np.nan)
        dat = sorted(unfold(dat_temp[i]))
        if i == 40:
            print("NNSD_temp = ",dat)
        # print(dat)
        for x in range(len_sub_list - 1):
            NN.append(dat[x + 1] - dat[x])
    # print(NN)
    return NN

def NNSD(data, mean_degree, hard, rate, eta=2):
    rad = str(mean_degree2rad(mean_degree, hard, rate, eta))
    list_len = len(data[rad])
    dat = [sorted(data[rad][i]) for i in range(list_len)]
    print("NNSD = ",dat[40])
    NN = [dat[i][x + 1] - dat[i][x] for i in range(list_len) for x in range(len(dat[i]) - 1)]
    # print(sorted(NN)[0:1000])
    return NN

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
dim = 1
rate = 1e3
mean_degree = 8

data = data_import(hard=hard, dim=dim, rate=rate, mean_degree=mean_degree)

NN = NNSD(data, mean_degree, hard, rate)
NN_temp = NNSD_temp(data, mean_degree, hard, rate)

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
