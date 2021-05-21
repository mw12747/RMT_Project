from RMTstats import random_geometric_graph, ray_conn_fun, hard_conn_fun, AdjSpectrum
import json
import os
from pathlib import Path
from scipy.special import gamma, gammainc
from math import sqrt, pi
from sys import platform


def run_sims_hard(rad, rate, dim, n_runs):

    n = 0
    eig_list = []
    print('rad = ', rad)

    while n < n_runs:
        if (n / n_runs * 100) % 10 == 0:
            print(n / n_runs * 100, "% complete")
        G, pos_nodes = random_geometric_graph(hard_conn_fun, rate=rate, dim=dim, torus=True, rad=rad)
        eig_list.append(AdjSpectrum(G))
        n += 1

    rad_key = str(rad)

    if platform == 'darwin':
        fname = f'/Users/michaelwilsher/Documents/GitHub/RMT_Project/SimulatedData/HardRGG/' \
                f'Dim_{dim}/rate_{rate}/eta_inf.txt'

    elif platform == 'win32':
        fname = f'C:/Users/mw12747/OneDrive - University of Bristol/MyFiles-Migrated/Documents/PycharmProjects/' \
                'RMT_Project/SimulatedData/HardRGG/Dim_{dim}/rate_{rate}/eta_inf.txt'

    path = fname[:-12]

    if not os.path.exists(path):
        os.makedirs(path)

    if os.path.isfile(fname):
        with open(fname, "r") as my_file:
            data = json.load(my_file)
        with open(fname, "w") as my_file:
            if rad_key in data:
                data[rad_key] += eig_list
            else:
                data[rad_key] = eig_list
            json.dump(data, my_file)

    else:
        with open(fname, "w") as my_file:
            data = {rad_key: eig_list}
            json.dump(data, my_file)


def run_sims_soft(rad, eta, rate, dim, n_runs):

    n = 0
    eig_list = []
    print('rad = ', rad)

    while n < n_runs:
        if (n / n_runs * 100) % 10 == 0:
            print(n / n_runs * 100, "% complete")
        H, pos_nodes = random_geometric_graph(ray_conn_fun, rate=rate, dim=dim, torus=True, rad=rad, eta=eta)
        eig_list.append(AdjSpectrum(H))
        n += 1

    rad_key = str(rad)

    if platform == 'darwin':
        fname = f'/Users/michaelwilsher/Documents/GitHub/RMT_Project/SimulatedData/SoftRGG/' \
                f'Dim_{dim}/rate_{rate}/eta_{eta}.txt'

    elif platform == 'win32':
        fname = f'C:/Users/mw12747/OneDrive - University of Bristol/MyFiles-Migrated/Documents/PycharmProjects/' \
                f'RMT_Project/SimulatedData/SoftRGG/Dim_{dim}/rate_{rate}/eta_{eta}.txt'

    path = fname[:-10]

    if not os.path.exists(path):
        os.makedirs(path)

    if os.path.isfile(fname):
        with open(fname, "r") as my_file:
            data = json.load(my_file)
        with open(fname, "w") as my_file:
            if rad_key in data:
                data[rad_key] += eig_list
            else:
                data[rad_key] = eig_list
            json.dump(data, my_file)

    else:
        with open(fname, "w") as my_file:
            data = {rad_key: eig_list}
            json.dump(data, my_file)


# def r_h(r_s, eta):
#     r_h = r_s / eta * (gamma(1 / eta) * gammainc(1 / eta, (2*r_s)**-eta))


def run_sims(n_runs, dim, mean_degree, rate):
    for d in mean_degree:
        print('mean_degree = ', d)
        if dim == 1:
            r_wax = d / (2 * rate)
            r_ray = d / (rate * sqrt(pi))
            r_hard = d / (2 * rate)
        elif dim == 2:
            r_wax = sqrt(d / (2 * pi * rate))
            r_ray = sqrt(d / (pi * rate))
            r_hard = sqrt(d / (pi * rate))

        run_sims_soft(r_wax, 1, rate, dim, n_runs)
        run_sims_soft(r_ray, 2, rate, dim, n_runs)
        run_sims_hard(r_hard, rate, dim, n_runs)


run_sims(n_runs=10, dim=1, mean_degree=[1], rate=1e3)

