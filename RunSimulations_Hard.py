from RMT_Project.RMTstats import random_geometric_graph, ray_conn_fun, hard_conn_fun, AdjSpectrum
import pandas as pd
import os
from pathlib import Path


def run_sims_hard(rad, rate, dim, n_runs):
    fname = '/Users/michaelwilsher/PycharmProjects/1DSoftRGG/RMT_Project/SimulatedData/HardRGG/Dim_%s/rate_%s/rad_%s/eta_inf.csv'\
            % (dim, rate, rad)

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

    fname_readme = '/Users/michaelwilsher/PycharmProjects/1DSoftRGG/RMT_Project/SimulatedData/HardRGG/' \
                   'Dim_%s/rate_%s/rad_%s/readme.txt' % (dim, rate, rad)

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


run_sims_hard(0.001, 1e3, 1, 100)
