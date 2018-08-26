from __future__ import division
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
import matplotlib.cm as cm
import matplotlib.colors as colors
from PIL import Image
from libtiff import TIFF
from osgeo import gdal
import ogr
from IPython import embed
import shapefile
from shapely.geometry import shape, Point
import matplotlib.patches as patches
from math import sin, cos, sqrt, atan2, radians, pi, degrees
from geopy.geocoders import Nominatim
geolocator = Nominatim()
import geopy.distance
from scipy import ndimage
from scipy.signal import convolve2d
from sklearn import linear_model
from scipy.interpolate import RectSphereBivariateSpline

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
    xOut=sqrt(xOut)
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

def scale(var):
	varScaled=(var-np.amin(var))/(np.amax(var)-np.amin(var))
	return varScaled

from itertools import compress
class MaskableList(list):
    def __getitem__(self, index):
        try: return super(MaskableList, self).__getitem__(index)
        except TypeError: return MaskableList(compress(self, index))

def checkIfInPolygon(lon, lat,polygon):
	point = Point(lon, lat)
	return polygon.contains(point)

def findGridDataLays(shapePoints,grid):
	'''Convert Vector To Raster'''
	griddedHits=np.zeros(shape=(len(grid[:,0]),len(grid[0,:])))
	for i in range(len(shapePoints)-1):
		x1,y1=shapePoints[i]
		x2,y2=shapePoints[i+1]
		if geopy.distance.distance([y1,x1],[y2,x2]).km<2:
			##### Find Grid Cell of Midpoint #####
			linePoints=np.array([(x1,y1),((x1+x2)/2.,(y1+y2)/2.),(x2,y2)])
		else:
			dist=geopy.distance.distance([y1,x1],[y2,x2]).km+1e-6
			if dist>50:
				continue
			if x2<x1:
				x=np.arange(x2,x1,abs(x2-x1)/(2*dist))
			else:
				if x2==x1:
					x2=x1+1e-4
				x=np.arange(x1,x2,abs(x2-x1)/(2*dist))
			slope,bInt=np.polyfit([x1,x2],[y1,y2],1)
			yfit=slope*np.array(x)+bInt
			linePoints=np.zeros(shape=(len(yfit),2))
			for j in range(len(yfit)):
				linePoints[j,:]=x[j],yfit[j]
		for j in range(len(linePoints)):
			midPoint=linePoints[j]
			yClosestL=min(grid[0,:,1], key=lambda i:abs(i-midPoint[1]))
			xClosestL=min(grid[:,0,0], key=lambda i:abs(i-midPoint[0]))
			if xClosestL<midPoint[0]:
				xClosestGrid=np.where(grid[:,0,0]==xClosestL)[0][0]
			else:
				xClosestGrid=np.where(grid[:,0,0]==xClosestL)[0][0]-1
			if yClosestL<midPoint[1]:
				yClosestGrid=np.where(grid[0,:,1]==yClosestL)[0][0]-1
			else:
				yClosestGrid=np.where(grid[0,:,1]==yClosestL)[0][0]
			griddedHits[xClosestGrid,yClosestGrid]=1
	return griddedHits # lon lat
				
		
def findGridAndDistance(shapePoints,grid):
	''' Give it a list of points in a road, will return pixel of center of each segment'''
	midpointHits=np.zeros(shape=(len(shapePoints)-1,2))
	dist=np.zeros(shape=(len(shapePoints)-1))
	for i in range(len(shapePoints)-1):
		x1,y1=shapePoints[i]
		x2,y2=shapePoints[i+1]
		##### Find Grid Cell of Midpoint #####
		midPoint=(x1+x2)/2.,(y1+y2)/2.

		yClosestL=min(grid[0,:,1], key=lambda i:abs(i-midPoint[1]))
		xClosestL=min(grid[:,0,0], key=lambda i:abs(i-midPoint[0]))
		if xClosestL<midPoint[0]:
			xClosestGrid=np.where(grid[:,0,0]==xClosestL)[0][0]
		else:
			xClosestGrid=np.where(grid[:,0,0]==xClosestL)[0][0]-1
		if yClosestL<midPoint[1]:
			yClosestGrid=np.where(grid[0,:,1]==yClosestL)[0][0]-1
		else:
			yClosestGrid=np.where(grid[0,:,1]==yClosestL)[0][0]
			midpointHits[i]=xClosestGrid,yClosestGrid

		### Convert radians to meters ###
		R = 6373.0 # Earth's radius
		coord1=y1,x1
		coord2=y2,x2
		dist[i]=geopy.distance.distance(coord1,coord2).km
		
	return midpointHits,dist

def createMidpointGrid(grid,pixelsize):
	gridMid=np.zeros(shape=(grid.shape))
	for ilon in range(len(grid[:,0,0])):
		gridMid[ilon,:,0]=np.mean([grid[ilon,0,0],grid[ilon,0,0]+pixelsize])
	for ilat in range(len(grid[0,:,1])):
		gridMid[:,ilat,1]=np.mean([grid[:,ilat,1],grid[:,ilat,1]+pixelsize])
	return gridMid

def findNearest(cityLonLat,gridMid,imageMask1):
	'''Give it:
	1. one dimentional lon lat data set
	2. grid of the mid points of cells (aquire from createMidpointGrid)
	3. mask for the grid in lon lat
	Will return:
	1. an array of the closest distance to one of the original lat lons (cityLonLat)
	2. which index of cityLonLat is closest
	'''
	lenLon,lenLat=len(gridMid[:,0,0]),len(gridMid[0,:,0])
	closestDist=10000*np.ones(shape=(lenLon,lenLat))
	iclosestDist=np.zeros(shape=(lenLon,lenLat),dtype=int)
	R = 6373.0 #earth's radius
	for ilon in range(lenLon):
		for ilat in range(lenLat):
			if imageMask1[ilon,ilat]==True:
				continue
			distances=np.sqrt((grid[ilon,ilat,0]-cityLonLat[:,0])**2+(grid[ilon,ilat,1]-cityLonLat[:,1])**2)
			leastDist=np.amin(distances)
			iclosestDist[ilon,ilat]=np.where(distances==leastDist)[0][0]
			closestDist[ilon,ilat]=geopy.distance.distance([gridMid[ilon,ilat,1],gridMid[ilon,ilat,0]],[cityLonLat[iclosestDist[ilon,ilat],1],cityLonLat[iclosestDist[ilon,ilat],0]]).km
		print np.round(100*ilon/float(lenLon),2),'%'

	closestDistM=np.ma.masked_array(closestDist,imageMask1)
	iclosestDistM=np.ma.masked_array(iclosestDist,imageMask1)
	closestDistM=np.swapaxes(closestDistM,0,1)
	iclosestDistM=np.swapaxes(iclosestDistM,0,1)
	return closestDistM,iclosestDistM

################################################
Lillian=True
Garyk=False

if Lillian:
	wddata='/Users/lilllianpetersen/iiasa/data/'
	wdfigs='/Users/lilllianpetersen/iiasa/figs/'
	wdvars='/Users/lilllianpetersen/iiasa/saved_vars/'

if Garyk:
	wddata='C:/Users/garyk/Documents/python_code/riskAssessmentFromPovertyEstimations/data/'
	wdfigs='C:/Users/garyk/Documents/python_code/riskAssessmentFromPovertyEstimations/figs/'
	wdvars='C:/Users/garyk/Documents/python_code/riskAssessmentFromPovertyEstimations/vars/'


makePlots=False

##############################################
# Gridded Malnutrition
##############################################
print 'Malnutrition'
tifWasting=TIFF.open(wddata+'malnutrition/IHME_AFRICA_CGF_2000_2015_WASTING_MEAN_2010_PREVALENCE_Y2018M02D28.TIF',mode='r')
ds=gdal.Open(wddata+'malnutrition/IHME_AFRICA_CGF_2000_2015_WASTING_MEAN_2015_PREVALENCE_Y2018M02D28.TIF')
width1 = ds.RasterXSize
height1 = ds.RasterYSize
gt = ds.GetGeoTransform()
minx = gt[0]
miny = gt[3] + width1*gt[4] + height1*gt[5] 
maxx = gt[0] + width1*gt[1] + height1*gt[2]
maxy = gt[3] 
pixelsize=abs(gt[-1])

# lat and lon
latm=np.ones(shape=(height1))
lonm=np.ones(shape=(width1))
for w in range(width1):
	lonm[w]=minx+w*pixelsize
for h in range(height1):
	latm[h]=miny+h*pixelsize
latm=latm[::-1] # reverse the order

mal=tifWasting.read_image()
mal[mal>1]=0.0308212 # The average of this weird pixel's neighbors

mal=mal[:1740]
latm=latm[:1740]
mal=mal[:,185:]
lonm=lonm[185:]

gridm=np.zeros(shape=(len(lonp),len(latp),2))
for x in range(len(lonm)):
	gridm[x,:,0]=lonm[x]
for y in range(len(latm)):
	gridm[:,y,1]=latm[y]

africaMask=np.array(mal)
africaMask[africaMask>=0]=0
africaMask[africaMask<0]=1

mal[mal<0]=0
plt.clf()
plt.imshow(mal,cmap=cm.jet,vmax=0.3)
plt.title('gridded wasting 2015')
plt.colorbar()
plt.savefig(wdfigs+'wasting',dpi=700)

##############################################
# Gridded Population
##############################################
print 'Population'
ds=gdal.Open(wddata+'population/gpw_v4_basic_demographic_characteristics_rev10_a000_004bt_2010_cntm_30_sec.tif')
width1 = ds.RasterXSize
height1 = ds.RasterYSize
gt = ds.GetGeoTransform()
minx = gt[0]
miny = gt[3] + width1*gt[4] + height1*gt[5] 
maxx = gt[0] + width1*gt[1] + height1*gt[2]
maxy = gt[3] 
pixelsize=abs(gt[-1])

# lat and lon
latp=np.ones(shape=(height1))
lonp=np.ones(shape=(width1))
for w in range(width1):
	lonp[w]=minx+w*pixelsize
for h in range(height1):
	latp[h]=miny+h*pixelsize
latp=latp[::-1] # reverse the order

worldPop=ds.ReadAsArray()
##### Scale to Africa #####
pop=worldPop[latp<np.amax(latm)+pixelsize]
latp=latp[latp<np.amax(latm)+pixelsize]
pop=pop[latp>np.amin(latm)]
latp=latp[latp>np.amin(latm)]

pop=pop[:,lonp<np.amax(lonm)+pixelsize]
lonp=lonp[lonp<np.amax(lonm)+pixelsize]
pop=pop[:,lonp>np.amin(lonm)]
lonp=lonp[lonp>np.amin(lonm)]

gridp=np.zeros(shape=(len(lonp),len(latp),2))
for x in range(len(lonp)):
	gridp[x,:,0]=lonp[x]
for y in range(len(latp)):
	gridp[:,y,1]=latp[y]

## Plot old and new grid
#plt.clf()
#plt.plot(np.ma.compressed(gridp[-8:,-8:,0]),np.ma.compressed(gridp[-8:,-8:,1]),'*k')
#plt.plot(np.ma.compressed(grid[-8:,-8:,0]),np.ma.compressed(grid[-8:,-8:,1]),'*r') 
#plt.savefig(wdfigs+'old_new_grid.pdf')
##

# This is the important line	
pop[pop<0]=0

plt.clf()
plt.imshow(pop,cmap=cm.gist_ncar_r,vmax=5000)
plt.title('gridded population')
plt.colorbar()
plt.savefig(wdfigs+'pop',dpi=700)

latl=np.radians(latp[::-1])+1.2
lonl=np.radians(lonp)+1.2
lutpop=RectSphereBivariateSpline(latl, lonl, pop)

newLats,newLons=np.meshgrid(np.radians(latm[::-1])+1.2,np.radians(lonm)+1.2)
pop1=lutpop.ev(newLats.ravel(),newLons.ravel()).reshape((len(lonm),len(latm))).T
pop1[pop1<0]=0

R=6371 #km
latdiff1=abs(np.sin(np.radians(latm[1:]))-np.sin(np.radians(latm[:-1])))
latdiff=np.zeros(shape=(len(latm)))
latdiff[:len(latdiff1)]=latdiff1
latdiff[-1]=latdiff[-2]

londiff1=abs(np.radians(lonm[1:])-np.radians(lonm[:-1]))
londiff=np.zeros(shape=(len(lonm)))
londiff[:len(londiff1)]=londiff1
londiff[-1]=londiff[-2]-(londiff[-3]-londiff[-2])

area=np.zeros(shape=(pop1.shape))
for ilon in range(len(londiff)):
	area[:,ilon]=(R**2)*latdiff*londiff[ilon]  #radians

pop5km=pop1*area
pop5km=np.ma.masked_array(pop5km,africaMask)

plt.clf()
plt.imshow(pop5km,cmap=cm.gist_ncar_r,vmax=1000)
plt.title('gridded population')
plt.colorbar()
plt.savefig(wdfigs+'pop5km',dpi=700)

malnumber=mal*pop5km
malnumber=np.ma.masked_array(malnumber,africaMask)

plt.clf()
plt.imshow(malnumber,cmap=cm.nipy_spectral_r,vmax=500)
plt.title('malnutrition number')
plt.colorbar()
plt.savefig(wdfigs+'malnumber',dpi=700)
exit()

######################################################
# Travel Time
######################################################
print 'Roads'

#rRoads = shapefile.Reader(wddata+'openstreetmap/'+country+'/openstreetmap/nga_trs_roads_osm.shp')
traveltif=TIFF.open(wddata+'travel_time/accessibility_to_cities_2015_v1.0.tif',mode='r')

ds=gdal.Open(wddata+'travel_time/accessibility_to_cities_2015_v1.0.tif')
width = ds.RasterXSize
height = ds.RasterYSize
gt = ds.GetGeoTransform()
minx = gt[0]
miny = gt[3] + width*gt[4] + height*gt[5] 
maxx = gt[0] + width*gt[1] + height*gt[2]
maxy = gt[3] 
pixelsizel=abs(gt[-1])

latt=np.ones(shape=(height))
lont=np.ones(shape=(width))
for w in range(width):
	lont[w]=minx+w*pixelsizel
for h in range(height):
	latt[h]=miny+h*pixelsizel

latt=latt[::-1]
if latm[0]<latm[-1]:
	latm=latm[::-1]

#travelWorld=traveltif.read_image()
##### Scale to Africa #####
travel=travelWorld[latt<np.amax(latm)+pixelsize]
latt=latt[latt<np.amax(latm)+pixelsize]
travel=travel[latt>np.amin(latm)]
latt=latt[latt>np.amin(latm)]

travel=travel[:,lont<np.amax(lonm)+pixelsize]
lont=lont[lont<np.amax(lonm)+pixelsize]
travel=travel[:,lont>np.amin(lonm)]
lont=lont[lont>np.amin(lonm)]

travel[travel<0]=-10

plt.clf()
plt.imshow(travel,cmap=cm.gist_ncar_r,norm=colors.LogNorm())
plt.title('Travel Time')
plt.colorbar()
plt.savefig(wdfigs+'travel',dpi=700)

### Cities ###

city=shapefile.Reader(wddata+'grump-v1-settlement-points-rev01-shp/global_settlement_points_v1.01.shp')

records=city.records()
citylonlat=np.zeros(shape=(5000,2))
cityname=[]
citycountry=[]
i=0
for record in records:
    if record[29]=='Africa':
        citylonlat[i,0]=record[3]
        citylonlat[i,1]=record[4]
        cityname.append(record[10])
        citycountry.append(record[30])
        i+=1




