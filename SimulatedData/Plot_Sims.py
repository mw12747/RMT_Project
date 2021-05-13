import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def plot_sims(hard, dim, rate, rad, eta=2):
    if hard:
        fname = '/Users/michaelwilsher/PycharmProjects/1DSoftRGG/' \
            'RMT_Project/SimulatedData/HardRGG/Dim_%s/rate_%s/rad_%s/eta_inf.csv' % (dim, rate, rad)
    else:
        fname = '/Users/michaelwilsher/PycharmProjects/1DSoftRGG/' \
                'RMT_Project/SimulatedData/SoftRGG/Dim_%s/rate_%s/rad_%s/eta_%s.csv' % (dim, rate, rad, eta)

    df = pd.read_csv(fname, index_col=False)

    data = {}
    for key in df.keys():
        data[key] = df[key].tolist()

    bins = np.arange(-40.05, 40.05, 0.1)

    for key in [str(rad)]:
        eig_temp = data[key]
        eig = [x for x in eig_temp if np.isnan(x) == False]
        print(len(eig))
        plt.hist(eig, bins=bins, density=True, label='rad = %s' % key)
        if hard:
            plt.title('Hard RGG for d = %s, rate = %s' % (dim, rate))
        else:
            plt.title('Soft RGG for d = %s, rate = %s, $\eta$ = %s' % (dim, rate, eta))


    plt.legend()
    plt.show()


plot_sims(True, 2, 1e3, 0.1, 1)
