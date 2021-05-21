# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 10:29:46 2015

New Edition Created on Monday Mar 22 15:05:00 2021

@author: Georgie & Michael
"""



"""
Here we calculate spectral statistics of RGGs
"""


"""
Import the necessary packages
"""

import random
import networkx as nx
import numpy as np
from scipy.spatial import distance
import scipy as sp



"""
Some formulas we will use 
"""


def random_geometric_graph(conn_fun, rate, dim, torus=True, **kwargs):

    def PP_nodes(rate, dim, *kargs):
        """
        Can be used at a later date to define a more general point process
        :param rate:
        :param dim:
        :param kargs:
        :return:
        """


    def PPP_nodes(rate, dim):
        """
        Returns a networkx Graph object with Poi(n_nodes) distributed uniformly at random in a dim-dimensional unit cube
        I.e. a dim-dimensional n_nodes rate PPP in a dim-dimensional unit cube
        :param rate: rate of the PPP used to generate the nodes in the RGG
        :param dim: [int] dimension of the system
        :return:
            n_nodes = [int/float] number of nodes in the node set
            pos_nodes = [list] positions of the nodes
        """

        # Check input types
        if type(dim) != int:
            raise TypeError("input parameter 'dim' must be an integer")

        n_nodes = np.random.poisson(rate)
        if dim == 1:
            pos_nodes = [[random.uniform(0, 1), 0] for n in range(n_nodes)]
        else:
            pos_nodes = [[random.random() for i in range(0,dim)] for n in range(n_nodes)]
        return n_nodes, pos_nodes


    def edges(rate, dim, torus=True):
        """
        Generates the edge set for a set of nodes distributed in dim-dimensional space according to the connection
        function, conn_fun
        :param dim: [int]
        :param torus: [bool]
        :return:
        """

        # Check input types
        if type(dim) != int:
            raise TypeError("input parameter 'dim' must be an integer")
        if type(torus) != bool:
            raise TypeError("input parameter 'torus' must be a boolean")

        n_nodes, pos_nodes = PPP_nodes(rate, dim)
        node_list = [i for i in range(n_nodes)]
        edge_list = []

        for i in range(n_nodes):
            for j in range(i + 1, n_nodes, 1):
                if torus:
                    # Need to generate toroidal distance
                    dist_temp = [min(abs(pos_nodes[i][d] - pos_nodes[j][d]), 1 - abs(pos_nodes[i][d] - pos_nodes[j][d]))
                                 for d in range(dim)]
                    dist = np.linalg.norm(dist_temp)
                else:
                    dist = distance.euclidean(pos_nodes[i], pos_nodes[j])
                if conn_fun(dist, **kwargs):
                        edge_list.append((i, j))

        return edge_list, node_list, pos_nodes

    edge_list, node_list, pos_nodes = edges(rate, dim, torus)

    G = nx.Graph()  # Create graph object G
    G.add_edges_from(edge_list)
    G.add_nodes_from(node_list)
    for n in range(len(node_list)):
        G.nodes[n]['pos'] = pos_nodes[n]

    return G, pos_nodes


def ray_conn_fun(dist, rad, eta):
    """
    Defines the Rayleigh connection function
    :param dist: distance between the nodes
    :param mu: inverse of the "connection range"
    :return: Returns the probability of nodes a distance dist apart being connected
    """
    prob = random.uniform(0, 1) < np.exp(-(dist / rad) ** eta)
    return prob


def hard_conn_fun(dist, rad):
    """
    Defines
    :param dist:
    :return:
    """
    prob = dist < rad
    return prob


def AdjSpectrum(G):
    """
    Calculate the eigenvalues of the adjacency matrix of a network
    :param G: [networkx graph object] this is a networkx graph object
    :return: [list] Returns a list of the real part of the eigenvalues of the adjacency matrix of G
    """
    eig_temp = nx.adjacency_spectrum(G)
    eig = [x.real for x in eig_temp]
    return eig


def LapSpectrum(G):
    """
    Calculate the eigenvalues of the Laplacian matrix of a network
    :param G: [networkx graph object] this is a networkx graph object
    :return: [list] Returns a list of the real part of the eigenvalues of the Laplacian matrix of G
    """
    eig_temp = nx.laplacian_spectrum(G)
    eig = [x.real for x in eig_temp]
    return eig


def Bincentres(arr):
    """ Return the mid points in an array """
    c = []
    for i in range(len(arr) - 1):
        c.append((arr[i] + arr[i + 1]) * 0.5)
    return c


def Spacings(ray):
    """Calculate the spacings between the values of an
    ordered array"""
    s = []
    for t in range(1, len(ray)):
        s.append(ray[t] - ray[t - 1])
    return s


def Spacings2(ray):
    """Calculate the next nearest neighbour spacings
    between the values of an ordered array"""
    s = []
    for t in range(len(ray) - 2):
        s.append((ray[t + 2] - ray[t]) * 0.5)
    return s


def Delta(vals, start, L):
    """Calculate the spectral rigidity for a sequence"""
    recentred = [i - (start + (L / 2.0)) for i in vals]
    n = len(recentred)
    recentred_sq = [i ** 2 for i in recentred]
    arr3 = []
    for j in range(n):
        arr3.append((n - (2.0 * (j + 1.0)) + 1.0) * (recentred[j]))
    return ((n ** 2) / (16.0)) - (1.0 / (L ** 2)) * (sum(recentred) ** 2) + ((3.0 * n) / (2.0 * (L ** 2))) * (
        sum(recentred_sq)) - (((3.0) / (L ** 4)) * (sum(recentred_sq) ** 2)) + (((1.0) / (L)) * (sum(arr3)))


def Brody(x, b):
    """ Returns the Brody distribution value of x with parameter b"""
    a = (sp.special.gamma((b + 2.0) / (b + 1.0))) ** (b + 1.0)
    p = (b + 1.0) * (a) * (x ** b) * (np.exp(-a * (x ** (b + 1.0))))
    return p


def ecdf(data):
    """Compute ECDF for a one-dimensional array of measurements."""
    # Number of data points: n
    n = len(data)
    # x-data for the ECDF: x
    x = np.sort(data)
    # y-data for the ECDF: y
    y = 1000.0 * (np.arange(1, n + 1) / n)
    return x, y


def GOPrig(x, c):
    """ Return the asymptotic ln function of the
    spectral rigidity of the GOE"""
    return ((1.0 / np.pi) ** 2) * (np.log(2.0 * np.pi * x) - 1.906) + c / x


def GSEnext(x):
    """a function which gives GSE data for nNNSD"""
    return 11.59745712271145 * (x ** 4) * (np.exp(-2.2635369684180673 * (x ** 2)))

# rad = 0.1
# rate = 100
# dim = 1
# eta = 2
#
# G, pos_nodes_g = random_geometric_graph(hard_conn_fun, rate, dim, torus=True)
#
# eta = 2
#
# H, pos_nodes_h = random_geometric_graph(ray_conn_fun, rate, dim, torus=True)
#
# eta = 1
#
# I, pos_nodes_i = random_geometric_graph(ray_conn_fun, rate, dim, torus=True)
#
# eig_G = AdjSpectrum(G)
# lap_eig_g = LapSpectrum(G)
# eig_H = AdjSpectrum(H)
# lap_eig_h = LapSpectrum(H)
# eig_I = AdjSpectrum(I)
# lap_eig_i = LapSpectrum(I)
#
# bins = np.arange(-60, 400, 0.5)
#
# plt.hist(eig_G, bins=bins, label='Hard RGG')
# plt.hist(eig_H, bins=bins, label='Soft RGG (Rayleigh)')
# plt.hist(eig_I, bins=bins, label='Soft RGG (Waxman)')
#
# plt.legend()
#
# plt.show()
#
# bins = np.arange(-10, 1000, 0.5)
#
# plt.hist(lap_eig_g, bins=bins, label='Hard RGG')
# plt.hist(lap_eig_h, bins=bins, label='Soft RGG (Rayleigh)')
# plt.hist(lap_eig_i, bins=bins, label='Soft RGG (Waxman)')
#
# plt.legend()
#
# plt.show()
