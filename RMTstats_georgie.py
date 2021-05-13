# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 10:29:46 2015

@author: Georgie
"""



"""
Here we calculate spectral statistics of RGGs
"""


"""
Import the necessary packages
"""

import math, random
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import time
import platform
import scipy as sp
import scipy.stats as stats
import scipy.optimize as opt
#import seaborn as sns
from scipy.optimize import curve_fit
from scipy import interpolate
#sns.set()




"""
Extract the name of the computer node. 
Set a clock to time the programme
"""

run=platform.node().split('.')[0]
t0=time.clock()/60.0

"""
Some formulas we will use 
"""


def RGG_nodes(nodes,dim): 
    """ returns a set of nodes randomly distributed in 
    the unit cube"""
    Nod=nx.Graph()
    Nod.add_nodes_from(range(nodes))
    for n in Nod:
        Nod.node[n]['pos']=[random.random() for i in range(0,dim)]
    return Nod    

def RGG_edges_Torus(Nodeset,r): 
    """take a set of nodes (Nodeset) and add edges to 
    create a RGG with connection radius r """
    nodes = Nodeset.nodes(data=True)
    while nodes:
        u,du = nodes.pop()
        pu = du['pos']
        for v,dv in nodes:
            pv = dv['pos']
            d = sum((min(abs(a-b),1.0-abs(a-b))**2 for a,b in zip(pu,pv)))
            if d <= r**2:
                Nodeset.add_edge(u,v)
    return Nodeset    


def SpecfinderAdj(Network):
    """ Calculate the adjacency spectrum of a network """
    a=nx.adjacency_spectrum(Network)
    ar=[t.real for t in a]
    return ar


def Bincentres(arr):
    """ Return the mid points in an array """
    c=[]
    for i in range(len(arr)-1):
        c.append((arr[i]+arr[i+1])*0.5)
    return c    
 
def Spacings(ray):
    """Calculate the spacings between the values of an 
    ordered array"""
    s=[]
    for t in range(1,len(ray)):
        s.append(ray[t]-ray[t-1])
    return s    

def Spacings2(ray):
    """Calculate the next nearest neighbour spacings 
    between the values of an ordered array"""
    s=[]
    for t in range(len(ray)-2):
        s.append((ray[t+2]-ray[t])*0.5)
    return s 
    
def Delta(vals,start,L):
    """Calculate the spectral rigidity for a sequence"""
    recentred = [i-(start+(L/2.0)) for i in vals]
    n=len(recentred)
    recentred_sq = [i**2 for i in recentred]
    arr3=[]
    for j in range(n):
        arr3.append((n-(2.0*(j+1.0))+1.0)*(recentred[j]))
    return ((n**2)/(16.0))-(1.0/(L**2))*(sum(recentred)**2) + ((3.0*n)/(2.0*(L**2)))*(sum(recentred_sq)) - (((3.0)/(L**4))*(sum(recentred_sq)**2)) + (((1.0)/(L))*(sum(arr3)))


def Brody(x,b):
    """ Returns the Brody distribution value of x with parameter b"""    
    a=(sp.special.gamma((b+2.0)/(b+1.0)))**(b+1.0)
    p=(b+1.0)*(a)*(x**b)*(np.exp(-a*(x**(b+1.0))))
    return p
    
def ecdf(data):
    """Compute ECDF for a one-dimensional array of measurements."""
    # Number of data points: n
    n = len(data)
    # x-data for the ECDF: x
    x = np.sort(data)
    # y-data for the ECDF: y
    y = 1000.0 * (np.arange(1, n+1)/n)
    return x, y
    
def GOPrig(x,c):
    """ Return the asymptotic ln function of the 
    spectral rigidity of the GOE"""
    return ((1.0/np.pi)**2)*(np.log(2.0*np.pi*x) -1.906)+c/x
    
def GSEnext(x):
    """a function which gives GSE data for nNNSD"""
    return 11.59745712271145*(x**4)*(np.exp(-2.2635369684180673*(x**2)))

   
"""
##################################################################
Parameters
"""

nod = 1000 #number of nodes
d = 2      #dimension
ra = 0.055   #connection radius
ens = 10000 #number of networks in the ensemble
stats_ens = 10000 #number of networks used to calculate the stats



"""
##################################################################
Create an ensemble of RGGs and save the ensembled spectrum to file.
Or load it from file (line 173)
"""

"""
aspec = [] # a list to store the adjacency spectrum

for nw in range(ens): # create an enseble and for each one extract the adjacency spectrum
    G=RGG_nodes(nod,d)
    RGG_edges_Torus(G,ra) 
    a=nx.adjacency_spectrum(G)
    ar=[t.real for t in a]
    aspec=aspec+ar
    #print(nw)
    
np.savetxt('AspecNod'+str(nod)+'Net'+str(ens)+'D'+str(d)+'R'+str(ra).replace(".","")+str(run)+'.txt',aspec)
"""

aspec = np.loadtxt('AspecNod1000Net10000D2R0055.txt')



"""
################################################################
Create a cumulative, normed histogram of the spectrum for use in 
unfolding
"""



sort_aspec = np.sort(aspec)
unfold = sp.interpolate.interp1d(sort_aspec,np.linspace(1.0, nod, sort_aspec.size),bounds_error=False, fill_value=np.nan)

#print("unfolder created")



"""
###################################################################
Create an ensemble of RGGs and calculate the statistics NNSD, nNNSD and Delta3
"""


NNSD = np.array([]) # array to store the nearest spacings
nNNSD = np.array([]) # array to store next nearest spacings
L_space = np.linspace(2.0,200,23)  # array of L values for Delta3

#D3_1 = []
#for j in range(23):
#    D3_1.append([])

#D3_2 = []
#for j in range(23):
#    D3_2.append([])    

D3 = []
for j in range(23):
    D3.append([])

for nwe in range(stats_ens):
    #print("Network # ", nwe)
    G = RGG_nodes(nod,d) # create a RGG
    RGG_edges_Torus(G,ra) 
    a = nx.adjacency_spectrum(G) # get its adjacency spectrum
    ar = [t.real for t in a]
    u = np.sort(unfold(ar)) # unfold it
    NNSD = np.append(NNSD, Spacings(u)) # Calculate the NNSD
    nNNSD = np.append(nNNSD, Spacings2(u)) # Calculate the nNNSD
    for ind, L in enumerate(L_space):
        start1 = np.arange(10,round(unfold([-1])[0]-L-10),L)
        start2 = np.arange(round(unfold([1])[0]+10),990-L,L)
        for s in start1:
            section = [eigenv for eigenv in u if  s <= eigenv < s+L]
            D3[ind].append(Delta(section, s, L))
        for s in start2:
            section = [eigenv for eigenv in u if  s <= eigenv < s+L]
            D3[ind].append(Delta(section, s, L))    

NNSD = NNSD[np.logical_not(np.isnan(NNSD))]
nNNSD = nNNSD[np.logical_not(np.isnan(nNNSD))]
D3_vals = [np.mean(v) for v in D3]

np.savetxt('NNSDNod'+str(nod)+'Net'+str(stats_ens)+'D'+str(d)+'R'+str(ra).replace(".","")+'.txt', NNSD)
np.savetxt('nNNSDNod'+str(nod)+'Net'+str(stats_ens)+'D'+str(d)+'R'+str(ra).replace(".","")+'.txt', nNNSD)
np.savetxt('D3Nod'+str(nod)+'Net'+str(stats_ens)+'D'+str(d)+'R'+str(ra).replace(".","")+'.txt', D3_vals)




"""
#############################################################
Load the data and plot
"""


"""
L_space = np.linspace(2.0,200,23)
#radii = ['003', '004','005', '006', '007', '008', '009', '01', '015', '02', '025', '03', '035', '04']
radii = ['003', '004', '005', '006', '04']
NNSD = {}
#nNNSD = {}
#D3 = {}
NNSDhist = {}
#BrodyValues = {}



for radius in radii:
    print(radius)
    NNSD[radius] = np.loadtxt('NNSDNod1000Net10000D2R'+radius+'.txt')
    #nNNSD[radius] = np.loadtxt('nNNSDNod1000Net10000D2R'+radius+'.txt')     
    #D3[radius] = np.loadtxt('D3Nod1000Net1000D2R'+radius+'.txt')
"""



"""
###################################################################
Drawing images for paper
"""






"""
#################################################################
NNSD
"""

"""
xvals = np.linspace(0.0,3.0, 10000)
GOE = [Brody(i,1.0) for i in xvals]
PO =  [Brody(i,0.0) for i in xvals]
plt.plot(xvals,GOE, color = 'black', linewidth = 1.0, label = 'GOE')
plt.plot(xvals,PO, '--', color = 'black', linewidth = 1.0, label = 'Poisson')


n1 = np.histogram(NNSD['003'][np.logical_and(NNSD['003'] > 0,  NNSD['003'] < 3.0)], bins = 100, normed = True)
xdata = Bincentres(n1[1])
ydata = n1[0]
plt.plot(xdata, ydata, 'o', label='r = 0.03', color = 'orange', alpha = 0.5)
popt, pcov = curve_fit(Brody, xdata, ydata,p0=(1.0))
br = [Brody(i,popt[0]) for i in xdata]
plt.plot(xdata, br, color = 'orange', linewidth = 0.8, alpha = 0.5)
print('Beta for 0.03 = ', popt[0])

n2 = np.histogram(NNSD['004'][np.logical_and(NNSD['004'] > 0,  NNSD['004'] < 3.0)], bins = 100, normed = True)
xdata = Bincentres(n2[1])
ydata = n2[0]
plt.plot(xdata, ydata, '<', label='r = 0.04', color = 'purple', alpha = 0.5)
popt, pcov = curve_fit(Brody, xdata, ydata,p0=(1.0))
br = [Brody(i,popt[0]) for i in xdata]
plt.plot(xdata, br, color = 'purple', linewidth = 0.8, alpha = 0.5)
print('Beta for 0.04 = ', popt[0])

n3 = np.histogram(NNSD['005'][np.logical_and(NNSD['005'] > 0,  NNSD['005'] < 3.0)], bins = 100, normed = True)
xdata = Bincentres(n3[1])
ydata = n3[0]
plt.plot(xdata, ydata, '.', label='r = 0.05', color = 'red', alpha = 0.5)
popt, pcov = curve_fit(Brody, xdata, ydata,p0=(1.0))
br = [Brody(i,popt[0]) for i in xdata]
plt.plot(xdata, br, color = 'red', linewidth = 0.8, alpha = 0.5)
print('Beta for 0.05 = ', popt[0])

n4 = np.histogram(NNSD['006'][np.logical_and(NNSD['006'] > 0,  NNSD['006'] < 3.0)], bins = 100, normed = True)
xdata = Bincentres(n4[1])
ydata = n4[0]
plt.plot(xdata, ydata, '>', label='r = 0.06', color = 'green', alpha = 0.5)
popt, pcov = curve_fit(Brody, xdata, ydata,p0=(1.0))
br = [Brody(i,popt[0]) for i in xdata]
plt.plot(xdata, br, color = 'green', linewidth = 0.8, alpha = 0.5)
print('Beta for 0.06 = ', popt[0])

n5 = np.histogram(NNSD['04'][np.logical_and(NNSD['04'] > 0,  NNSD['04'] < 3.0)], bins = 100, normed = True)
xdata = Bincentres(n5[1])
ydata = n5[0]
plt.plot(xdata, ydata, '*', label='r = 0.4', color = 'blue', alpha = 0.5)
popt, pcov = curve_fit(Brody, xdata, ydata,p0=(1.0))
br = [Brody(i,popt[0]) for i in xdata]
plt.plot(xdata, br, color = 'blue', linewidth = 0.8, alpha = 0.5)
print('Beta for 0.4 = ', popt[0])

plt.legend(fontsize = 20)
plt.xlim(0,3)
plt.xlabel('s',fontsize=20)
plt.ylabel('P(s)',fontsize=20)
txt='(a)'
plt.text(2.5,0.2,txt,fontsize=20)
plt.tight_layout()




#radius = radii[0]
#hist = np.histogram(NNSD[radius][NNSD[radius]>0], bins ='auto', normed = True)
#xdata = Bincentres(hist[1])
#ydata = hist[0]
#plt.plot(xdata, ydata, '.')
#popt, pcov = curve_fit(Brody, xdata, ydata,p0=(1.0))
#BrodyValues[radius] = popt
#br = [Brody(i,popt[0]) for i in xdata]
#plt.plot(xdata, br)





"""

"""
################################################################
Brody fit
"""

"""
for radius in radii:
    NNSDhist[radius] = np.histogram(NNSD[radius][np.logical_and(NNSD[radius] > 0,  NNSD[radius] < 3.0)], bins = 100 , normed = True)
    xdata = Bincentres(NNSDhist[radius][1])
    ydata = NNSDhist[radius][0]
    popt, pcov = curve_fit(Brody, xdata, ydata,p0=(1.0))
    BrodyValues[float(radius[0]+'.'+radius[1:])] = popt[0]


plt.plot(list(BrodyValues.keys()),list(BrodyValues.values()), 'o')
plt.xlim(0.025,0.41)
plt.ylim(0.0,1.0)
plt.xlabel('r',fontsize=20)
plt.ylabel(r'$\beta$',fontsize=20)
txt='(b)'
plt.text(0.35,0.7,txt,fontsize=20)
plt.tight_layout()
"""


"""
###############################################################
nNNSD
"""

"""
xvals = np.linspace(0.0,2.0, 10000)
GOE = [GSEnext(i) for i in xvals]
#PO =  [Brody(i,0.0) for i in xvals]
plt.plot(xvals,GOE, color = 'black', linewidth = 1.0, alpha = 0.6, label = 'GOE')
#plt.plot(xvals,PO, '--', color = 'black', linewidth = 1.0, label = 'Poisson')

n1 = np.histogram(nNNSD['003'][np.logical_and(nNNSD['003'] > 0,  nNNSD['003'] < 2.0)], bins = 75, normed = True)
xdata = Bincentres(n1[1])
ydata = n1[0]
plt.plot(xdata, ydata, 'o', label='r = 0.03', color = 'orange', alpha = 0.5)


n2 = np.histogram(nNNSD['004'][np.logical_and(nNNSD['004'] > 0,  nNNSD['004'] < 2.0)], bins = 75, normed = True)
xdata = Bincentres(n2[1])
ydata = n2[0]
plt.plot(xdata, ydata, '<', label='r = 0.04', color = 'purple', alpha = 0.5)

n3 = np.histogram(nNNSD['005'][np.logical_and(nNNSD['005'] > 0,  nNNSD['005'] < 2.0)], bins = 75, normed = True)
xdata = Bincentres(n3[1])
ydata = n3[0]
plt.plot(xdata, ydata, '.', label='r = 0.05', color = 'red', alpha = 0.5)


n4 = np.histogram(nNNSD['006'][np.logical_and(nNNSD['006'] > 0,  nNNSD['006'] < 2.0)], bins = 75, normed = True)
xdata = Bincentres(n4[1])
ydata = n4[0]
plt.plot(xdata, ydata, '>', label='r = 0.06', color = 'green', alpha = 0.5)


n5 = np.histogram(nNNSD['04'][np.logical_and(nNNSD['04'] > 0,  nNNSD['04'] < 2.0)], bins = 75, normed = True)
xdata = Bincentres(n5[1])
ydata = n5[0]
plt.plot(xdata, ydata, '*', label='r = 0.4', color = 'blue', alpha = 0.5)


plt.legend(fontsize = 20)
plt.xlim(0,2)
plt.xlabel(r'$s_2$',fontsize=20)
plt.ylabel(r'$P(s_2)$',fontsize=20)
#txt='(a)'
#plt.text(2.5,0.2,txt,fontsize=20)
plt.tight_layout()
"""



"""
##################################################################
Spectral density draw
"""

"""
s1 = np.loadtxt('AspecNod1000Net10000D2R01.txt')
s2 = np.loadtxt('AspecNod1000Net10000D2R03.txt')

plt.figure(1)
plt.hist(s1, bins = 'auto', normed = 'True',histtype = 'stepfilled')
plt.xlabel(r'$\lambda$', fontsize = 20)
plt.ylabel(r'$\rho(\lambda)$', fontsize = 20)
plt.xlim(-15,15)
txt='(a)'
plt.text(-10,0.8,txt,fontsize=20)
plt.tight_layout()


plt.figure(2)
plt.hist(s2, bins = 'auto', normed = 'True', histtype = 'stepfilled')
plt.xlabel(r'$\lambda$', fontsize = 20)
plt.ylabel(r'$\rho(\lambda)$', fontsize = 20)
plt.xlim(-40,40)
txt='(b)'
plt.text(-30,0.14,txt,fontsize=20)
plt.tight_layout()

"""






"""
###################################################################
Draw the NNSD difference to GOE and Brody difference to GOE
"""


"""
n1 = np.loadtxt('NNSDNod1000Net10000D2R005.txt') #NNSD 0.05
n2 = np.loadtxt('NNSDNod1000Net10000D2R04.txt') #NNSD 0.4
h1 = np.histogram(n1[np.logical_and(n1 > 0,  n1 < 3.0)], bins = 100 , normed = True)
h2 = np.histogram(n2[np.logical_and(n2 > 0,  n2 < 3.0)], bins = 100 , normed = True)
GOE1=[Brody(i,1.0) for i in Bincentres(h1[1])]
GOE2=[Brody(i,1.0) for i in Bincentres(h2[1])]
d1=GOE1-h1[0]
d2=GOE2-h2[0]

plt.figure(1)
plt.plot(Bincentres(h1[1]),d1,'.',color='red',label='r = 0.05')

plt.figure(2)
plt.plot(Bincentres(h2[1]),d2,'.',color='red',label='r = 0.4')


v=np.linspace(0.0,3.0,10000)
G1=np.array([Brody(i,1.0) for i in v])
B1=np.array([Brody(i,0.69561115) for i in v])
B2=np.array([Brody(i,0.95806816) for i in v])
bd1=G1-B1
bd2=G1-B2
plt.figure(1)
plt.plot(v,bd1,color='black')
plt.xlabel('$s$',fontsize=20)
plt.ylabel('$GOE-P(s)$',fontsize=20)
txt='(a)'
plt.text(2.5,0.1,txt,fontsize=20)
plt.tight_layout()

plt.figure(2)
plt.plot(v,bd2,color='black')
plt.xlabel('$s$',fontsize=20)
plt.ylabel('$GOE-P(s)$',fontsize=20)
txt='(b)'
plt.text(2.5,0.015,txt,fontsize=20)
plt.tight_layout()

"""








"""
###################################################################
Plot cumulative histogram
"""



"""
x,y = ecdf(aspec)

plt.step(x,y, color='blue')

G=RGG_nodes(nod,d)
RGG_edges_Torus(G,ra) 
a=nx.adjacency_spectrum(G)
ar=[t.real for t in a]
x2,y2 = ecdf(ar)

plt.step(x2,y2, color = 'red')
plt.xlabel('$\lambda$',fontsize=20)
plt.ylabel('$\eta(\lambda)$',fontsize=20)

a = plt.axes([.4, .2, .4, .4])
plt.step(x,y, color='blue')
plt.step(x2,y2, color = 'red')
plt.xlim(-5,-4)
plt.ylim(70,110)
plt.xticks([-5,-4])
plt.yticks([70,110])
"""


"""
#########################################################
SR
"""



"""
x = np.linspace(1.0,200,1000)
g = [GOPrig(i,0.05) for i in x]
p=[i/15.0 for i in x]
pf=[1/12.0 for i in x]
plt.plot(x,p,'--',color='green',label='Poisson')
plt.plot(x,pf,'-.',color='black',label='Picket Fence')
plt.plot(x,g,color='black',label='GOE')
plt.legend(loc = 7, fontsize = 15)
plt.plot(L_space, D3['005'], 'o', color = 'red', alpha = 0.7, label = 'r = 0.05')
plt.plot(L_space, D3['006'], 'd', color = 'orange', alpha = 0.7, label = 'r = 0.06')
plt.plot(L_space, D3['007'], 'D', color = 'blue', alpha = 0.7, label = 'r = 0.07')
plt.plot(L_space, D3['008'], '^', color = 'green',alpha = 0.7, label = 'r = 0.08')
plt.plot(L_space, D3['01'], 'v', color = 'red',alpha = 0.7, label = 'r=0.1')
plt.plot(L_space, D3['015'], 'p', color = 'orange',alpha = 0.7, label = 'r=0.15')
plt.plot(L_space, D3['02'], '.', color = 'blue',alpha = 0.7, label = 'r=0.2')
plt.plot(L_space, D3['04'], '*', color = 'green',alpha = 0.7, label = 'r = 0.4')
plt.ylim(0.0,0.8)
plt.xlabel('L',fontsize=20)
plt.ylabel('$\Delta_3(L)$',fontsize=20)
plt.tight_layout()
"""











"""
#########################################################
Draw images for poster
"""

"""
#NNSD
xvals = np.linspace(0.0,3.0, 10000)
GOE = [Brody(i,1.0) for i in xvals]
PO =  [Brody(i,0.0) for i in xvals]
plt.plot(xvals,GOE, color = 'black', linewidth = 1.0, label = 'GOE')
plt.plot(xvals,PO, '--', color = 'black', linewidth = 1.0, label = 'Poisson')

n1 = np.histogram(NNSD['005'][np.logical_and(NNSD['005'] > 0,  NNSD['005'] < 3.0)], bins = 'auto', normed = True)
xdata = Bincentres(n1[1])
ydata = n1[0]
plt.plot(xdata, ydata, '.', label='r=0.05', color = 'red', alpha = 0.5)
popt, pcov = curve_fit(Brody, xdata, ydata,p0=(1.0))
br = [Brody(i,popt[0]) for i in xdata]
plt.plot(xdata, br, color = 'red', linewidth = 0.8, alpha = 0.5)
print('Beta for 0.05 = ', popt[0])

n2 = np.histogram(NNSD['006'][np.logical_and(NNSD['006'] > 0,  NNSD['006'] < 3.0)], bins = 'auto', normed = True)
xdata = Bincentres(n2[1])
ydata = n2[0]
plt.plot(xdata, ydata, '.', label='r=0.06', color = 'green', alpha = 0.5)
popt, pcov = curve_fit(Brody, xdata, ydata,p0=(1.0))
br = [Brody(i,popt[0]) for i in xdata]
plt.plot(xdata, br, color = 'green', linewidth = 0.8, alpha = 0.5)
print('Beta for 0.06 = ', popt[0])

n3 = np.histogram(NNSD['025'][np.logical_and(NNSD['025'] > 0,  NNSD['025'] < 3.0)], bins ='auto', normed = True)
xdata = Bincentres(n3[1])
ydata = n3[0]
plt.plot(xdata, ydata, '.', label='r=0.25', color = 'blue', alpha = 0.5)
popt, pcov = curve_fit(Brody, xdata, ydata,p0=(1.0))
br = [Brody(i,popt[0]) for i in xdata]
plt.plot(xdata, br, color = 'blue', linewidth = 0.8, alpha = 0.5)
print('Beta for 0.25 = ', popt[0])

plt.legend()
plt.xlim(0,3)
plt.xlabel('s',fontsize=20)
plt.ylabel('$P(s)$',fontsize=20)

"""

# SR

"""
x = np.linspace(1.0,200,1000)
g = [GOPrig(i,0.05) for i in x]
p=[i/15.0 for i in x]
pf=[1/12.0 for i in x]
plt.plot(x,p,'--',color='green',label='Poisson')
plt.plot(x,pf,'-.',color='black',label='Picket Fence')
plt.plot(x,g,color='black',label='GOE')
plt.legend(loc = 7, fontsize = 15)
plt.plot(L_space, D3['005'], '.', color = 'red', alpha = 0.7, label = 'r = 0.05')
#plt.plot(L_space, D3['006'], '.', color = 'green',alpha = 0.7, label = 'r=0.06')
#plt.plot(L_space, D3['007'], '.', color = 'orange',alpha = 0.7, label = 'r=0.07')
plt.plot(L_space, D3['008'], '*', color = 'green',alpha = 0.7, label = 'r = 0.08')
#plt.plot(L_space, D3['009'], '.', color = 'brown',alpha = 0.7, label = 'r=0.09')
#plt.plot(L_space, D3['01'], '.', color = 'blue',alpha = 0.7, label = 'r = 0.1')
plt.plot(L_space, D3['015'], '>', color = 'blue',alpha = 0.7, label = 'r=0.15')
#plt.plot(L_space, D3['02'], '.', color = 'red',alpha = 0.7, label = 'r=0.2')
#plt.plot(L_space, D3['025'], '.', color = 'green',alpha = 0.7, label = 'r=0.25')
#plt.plot(L_space, D3['03'], '.', color = 'orange',alpha = 0.7, label = 'r=0.3')
#plt.plot(L_space, D3['035'], '.', color = 'grey',alpha = 0.7, label = 'r=0.35')
plt.plot(L_space, D3['04'], '<', color = 'orange',alpha = 0.7, label = 'r = 0.4')
plt.ylim(0.0,0.8)
plt.xlabel('L',fontsize=20)
plt.ylabel('$\Delta_3(L)$',fontsize=20)
plt.tight_layout()
"""

########################################
output=str(time.clock()/60.0-t0)+'\n'
print(output) 
########################################