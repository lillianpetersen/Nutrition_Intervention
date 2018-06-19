import csv
from math import sqrt
from sys import exit
import pickle
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.mlab as mlab
from mpl_toolkits.basemap import Basemap
import os
from scipy.stats import norm
import matplotlib as mpl
from matplotlib.patches import Polygon
import random

###############################################
# Functions
###############################################
def stdDev(x):
    '''function to compute standard deviation'''
    xAvg=np.mean(x)
    xOut=0.
    for k in range(len(x)):
        xOut=xOut+(x[k]-xAvg)**2
    xOut=xOut/(k+1)
    xOut=math.sqrt(xOut)
    return xOut

def Variance(x):
    '''function to compute the variance (std dev squared)'''
    xAvg=np.mean(x)
    xOut=0.
    for k in range(len(x)):
        xOut=xOut+(x[k]-xAvg)**2
    xOut=xOut/(k+1)
    return xOut

def SumOfSquares(x):
    '''function to compute the sum of squares'''
    xOut=0.
    for k in range(len(x)):
        xOut=xOut+x[k]**2
    return xOut

def corr(x,y):
    ''' function to find the correlation of two arrays'''
    xAvg=np.mean(x)
    Avgy=np.mean(y)
    rxy=0.
    n=min(len(x),len(y))
    for k in range(n):
        rxy=rxy+(x[k]-xAvg)*(y[k]-Avgy)
    rxy=rxy/(k+1)
    stdDevx=stdDev(x)
    stdDevy=stdDev(y)
    rxy=rxy/(stdDevx*stdDevy)
    return rxy

###############################################

wddata='/Users/lilllianpetersen/iiasa/data/'
wdfigs='/Users/lilllianpetersen/iiasa/figs/'
plt.clf()

################## Ethiopia OpenStreetMap ##################
# Ethiopia lat lons
map = Basemap(llcrnrlon=33,llcrnrlat=2.5,urcrnrlon=48.5,urcrnrlat=15.5, projection='lcc',lon_0=40.49,lat_0=9.145)

map.readshapefile(wddata+'ethiopia-openstreetmap/gis_osm_roads_free_1', name='roads', drawbounds=True)

ax = plt.gca()
plt.title('Ethiopian Roads from Open Street Map')
plt.yticks([])
plt.xticks([])

plt.savefig(wdfigs+'openstreetmap_ethiopia',dpi=700)
#poly = Polygon(seg,facecolor=[cmapArray[icmap,0],cmapArray[icmap,1],cmapArray[icmap,2]],edgecolor=[0,0,0])


################## Africa Rivers ##################

map = Basemap(llcrnrlon=33,llcrnrlat=2.5,urcrnrlon=48.5,urcrnrlat=15.5, projection='lcc',lon_0=40.49,lat_0=9.145)

map.readshapefile(wddata+'ethiopia-openstreetmap/gis_osm_roads_free_1', name='roads', drawbounds=True)

ax = plt.gca()
plt.title('Ethiopian Roads from Open Street Map')
plt.yticks([])
plt.xticks([])

plt.savefig(wdfigs+'openstreetmap_ethiopia',dpi=700)














