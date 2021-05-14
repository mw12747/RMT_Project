from RMTstats import random_geometric_graph, ray_conn_fun, hard_conn_fun, AdjSpectrum
import pandas as pd
import os
from pathlib import Path
from scipy.special import gamma, gammainc
from math import sqrt, pi
from sys import platform


def run_sims_hard(rad, rate, dim, n_runs):
    if platform == 'darwin':
        fname = '/Users/michaelwilsher/PycharmProjects/1DSoftRGG/RMT_Project/SimulatedData/HardRGG/Dim_%s/rate_%s/rad_%s/eta_inf.csv'\
                % (dim, rate, rad)

    if platform == 'win32':
        fname = 'C:/Users/mw12747/OneDrive - University of Bristol/MyFiles-Migrated/Documents/PycharmProjects/' \
                'RMT_Project/SimulatedData/HardRGG/Dim_%s/rate_%s/rad_%s/eta_inf.csv' % (dim, rate, rad)

    path = fname[:-12]

    if not os.path.exists(path):
        os.makedirs(path)

    if os.path.isfile(fname):
        df_old = pd.read_csv(fname, index_col=False)
        data = {}
        for key in df_old.keys():
            data[key] = df_old[key].tolist()
    else:
        data = {}

    if platform == 'darwin':
        fname_readme = '/Users/michaelwilsher/PycharmProjects/1DSoftRGG/RMT_Project/SimulatedData/HardRGG/' \
                       'Dim_%s/rate_%s/rad_%s/readme.txt' % (dim, rate, rad)

    if platform == 'win32':
        fname_readme = 'C:/Users/mw12747/OneDrive - University of Bristol/MyFiles-Migrated/Documents/PycharmProjects/' \
                       'RMT_Project/SimulatedData/HardRGG/Dim_%s/rate_%s/rad_%s/readme.txt' % (dim, rate, rad)

    if not os.path.isfile(fname_readme):
        Path(fname_readme).touch()
        file = open(fname_readme, "w")
        file.write("This folder contains the simulation data for a Hard RGG in Dimension = %s, Rate = %s, Radius = %s \n nsims = 0"
                   % (dim, rate, rad))
        nsims_old = 0
    else:
        file = open(fname_readme, "r")
        for line in file.readlines():
            if 'nsims' in line:
                nsims_temp = line[8:]
                nsims_old = float(nsims_temp)


    n = 0
    eig_list = []
    print('rad = ', rad)

    while n < n_runs:
        nsims_new = nsims_old + n
        if (n / n_runs * 100) % 10 == 0:
            print(n / n_runs * 100, "% complete")
        G, pos_nodes = random_geometric_graph(hard_conn_fun, rate=rate, dim=dim, torus=True, rad=rad)
        eig_list = eig_list + AdjSpectrum(G)
        n += 1

    rad_key = str(rad)

    if rad_key in data:
        data[rad_key] = data[rad_key] + eig_list
    else:
        data[rad_key] = eig_list

    df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in data.items()]))

    if not os.path.exists(path):
        os.makedirs(path)

    df.to_csv(fname, index=False)

    file = open(fname_readme, "w")
    file.write("This folder contains the simulation data for a Hard RGG in Dimension = %s, Rate = %s, Radius = %s \n nsims = %s"
               % (dim, rate, rad, nsims_new))


def run_sims_soft(rad, eta, rate, dim, n_runs):
    print(platform)
    if platform == 'darwin':
        fname = '/Users/michaelwilsher/PycharmProjects/1DSoftRGG/RMT_Project/SimulatedData/SoftRGG/Dim_%s/rate_%s/rad_%s/eta_%s.csv'\
            % (dim, rate, rad, eta)

    if platform == 'win32':
        fname = 'C:/Users/mw12747/OneDrive - University of Bristol/MyFiles-Migrated/Documents/PycharmProjects/' \
                'RMT_Project/SimulatedData/SoftRGG/Dim_%s/rate_%s/rad_%s/eta_%s.csv'\
            % (dim, rate, rad, eta)

    path = fname[:-10]

    if not os.path.exists(path):
        os.makedirs(path)

    if os.path.isfile(fname):
        df_old = pd.read_csv(fname, index_col=False)
        data = {}
        for key in df_old.keys():
            data[key] = df_old[key].tolist()
    else:
        data = {}

    if platform == 'darwin':
        fname_readme = '/Users/michaelwilsher/PycharmProjects/1DSoftRGG/RMT_Project/SimulatedData/SoftRGG/' \
                   'Dim_%s/rate_%s/rad_%s/readme_eta_%s.txt' % (dim, rate, rad, eta)

    if platform == 'win32':
        fname_readme = 'C:/Users/mw12747/OneDrive - University of Bristol/MyFiles-Migrated/Documents/PycharmProjects/' \
                'RMT_Project/SimulatedData/SoftRGG/' \
                   'Dim_%s/rate_%s/rad_%s/readme_eta_%s.txt' % (dim, rate, rad, eta)

    if not os.path.isfile(fname_readme):
        Path(fname_readme).touch()
        file = open(fname_readme, "w")
        file.write("This folder contains the simulation data for a Soft RGG in Dimension = %s, Rate = %s, Radius = %s, Eta = %s \n nsims = 0"
                   % (dim, rate, rad, eta))
        nsims_old = 0
    else:
        file = open(fname_readme, "r")
        for line in file.readlines():
            if 'nsims' in line:
                nsims_temp = line[8:]
                nsims_old = float(nsims_temp)


    n = 0
    eig_list = []
    print('rad = ', rad)

    while n < n_runs:
        nsims_new = nsims_old + n
        if (n / n_runs * 100) % 10 == 0:
            print(n / n_runs * 100, "% complete")
        H, pos_nodes = random_geometric_graph(ray_conn_fun, rate=rate, dim=dim, torus=True, rad=rad, eta=eta)
        eig_list = eig_list + AdjSpectrum(H)
        n += 1

    rad_key = str(rad)

    if rad_key in data:
        data[rad_key] = data[rad_key] + eig_list
    else:
        data[rad_key] = eig_list

    df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in data.items()]))

    if not os.path.exists(path):
        os.makedirs(path)

    df.to_csv(fname, index=False)

    file = open(fname_readme, "w")
    file.write("This folder contains the simulation data for a Soft RGG in Dimension = %s, Rate = %s, Radius = %s, Eta = %s \n nsims = %s"
               % (dim, rate, rad, eta, nsims_new))


def r_h(r_s, eta):
    r_h = r_s / eta * (gamma(1 / eta) * gammainc(1 / eta, (2*r_s)**-eta))


def run_sims(n_runs, dim, mean_degree, rate):
    for d in mean_degree:
        print('mean_degree = ', d)
        r_wax = d / (2 * rate)
        r_ray = d / (rate * sqrt(pi))
        r_hard = d / (2 * rate)
        run_sims_soft(r_wax, 1, rate, dim, n_runs)
        run_sims_soft(r_ray, 2, rate, dim, n_runs)
        run_sims_hard(r_hard, rate, dim, n_runs)

run_sims(n_runs=100, dim=1, mean_degree=[6,7,8,9,10], rate=1e3)

