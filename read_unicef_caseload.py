from __future__ import division
import csv
from math import sqrt
from sys import exit
import pickle
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.mlab as mlab
#from mpl_toolkits.basemap import Basemap
import os
from scipy.stats import norm, mode
import matplotlib as mpl
from matplotlib.patches import Polygon
import random
import matplotlib.cm as cm
import matplotlib.colors as colors
from PIL import Image
from osgeo import gdal
#import ogr
from IPython import embed
import shapefile
#from shapely.geometry import shape, Point
import matplotlib.patches as patches
from math import sin, cos, sqrt, atan2, radians, pi, degrees
from geopy.geocoders import Nominatim
geolocator = Nominatim()
import geopy.distance
from scipy import ndimage
from scipy.signal import convolve2d
from sklearn import linear_model
from scipy.interpolate import RectSphereBivariateSpline
from scipy.interpolate import griddata
from sklearn.ensemble import RandomForestRegressor
from scipy.interpolate import RegularGridInterpolator
import sklearn
import csv
#import googlemaps

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

def findGridDataLays(shapePoints,gridMid):
	'''Convert Vector To Raster'''
	griddedHits=np.zeros(shape=(len(gridMid[:,0]),len(gridMid[0,:])))
	stepsize=np.zeros(shape=(len(shapePoints)))
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
				step=abs(x2-x1)/(2*dist)
				if step==0:
					step==0.004
				x=np.arange(x2,x1,step)
			else:
				if x2==x1:
					x2=x1+1e-4
				step=abs(x2-x1)/(2*dist)
				if step==0:
					step==0.004
				x=np.arange(x1,x2,step)
			stepsize[i]=abs(x2-x1)/(2*dist)
			slope,bInt=np.polyfit([x1,x2],[y1,y2],1)
			yfit=slope*np.array(x)+bInt
			linePoints=np.zeros(shape=(len(yfit),2))
			for j in range(len(yfit)):
				linePoints[j,:]=x[j],yfit[j]
		for j in range(len(linePoints)):
			midPoint=linePoints[j]
			yClosestL=min(gridMid[0,:,1], key=lambda i:abs(i-midPoint[1]))
			xClosestL=min(gridMid[:,0,0], key=lambda i:abs(i-midPoint[0]))
			xClosestGrid=np.where(gridMid[:,0,0]==xClosestL)[0][0]
			yClosestGrid=np.where(gridMid[0,:,1]==yClosestL)[0][0]
			griddedHits[xClosestGrid,yClosestGrid]=1
	return griddedHits # lon lat
				
		
def findGridAndDistance(shapePoints,gridMid):
	''' Give it a list of points in a road, will return pixel of center of each segment'''
	midpointHits=np.zeros(shape=(len(shapePoints)-1,2))
	dist=np.zeros(shape=(len(shapePoints)-1))
	for i in range(len(shapePoints)-1):
		x1,y1=shapePoints[i]
		x2,y2=shapePoints[i+1]
		##### Find Grid Cell of Midpoint #####
		midPoint=(x1+x2)/2.,(y1+y2)/2.

		yClosestL=min(gridMid[0,:,1], key=lambda i:abs(i-midPoint[1]))
		xClosestL=min(gridMid[:,0,0], key=lambda i:abs(i-midPoint[0]))

		xClosestGrid=np.where(gridMid[:,0,0]==xClosestL)[0][0]
		yClosestGrid=np.where(gridMid[0,:,1]==yClosestL)[0][0]
		midpointHits[i]=xClosestGrid,yClosestGrid

		### Convert radians to meters ###
		R = 6373.0 # Earth's radius
		coord1=x1,y1
		coord2=x2,y2
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
			distances=np.sqrt((gridMid[ilon,ilat,0]-cityLonLat[:,0])**2+(gridMid[ilon,ilat,1]-cityLonLat[:,1])**2)
			leastDist=np.amin(distances)
			iclosestDist[ilon,ilat]=np.where(distances==leastDist)[0][0]
			closestDist[ilon,ilat]=geopy.distance.distance([gridMid[ilon,ilat,1],gridMid[ilon,ilat,0]],[cityLonLat[iclosestDist[ilon,ilat],1],cityLonLat[iclosestDist[ilon,ilat],0]]).km
		print np.round(100*ilon/float(lenLon),2),'%'

	closestDistM=np.ma.masked_array(closestDist,imageMask1)
	iclosestDistM=np.ma.masked_array(iclosestDist,imageMask1)
	return closestDistM,iclosestDistM

def marketPotentials(pop,popCoord,gridMid,imageMask2):
	lenLat=len(gridMid[:,0,0])
	lenLon=len(gridMid[0,:,0])
	MP=np.zeros(shape=(gridMid[:,:,0].shape))
	for ilat in range(lenLat):
		print np.round(100*ilat/float(lenLat),2),'%'
		for ilon in range(lenLon):
			if imageMask2[ilat,ilon]==True:
				continue
			dists=np.sqrt((gridMid[ilat,ilon,0]-popCoord[:,0])**2+(gridMid[ilat,ilon,1]-popCoord[:,1])**2)
			MP[ilat,ilon]=np.sum(pop/(dists**1.2))
	return MP

################################################
try:
	wddata='/Users/lilllianpetersen/iiasa/data/'
	wdfigs='/Users/lilllianpetersen/iiasa/figs/'
	wdvars='/Users/lilllianpetersen/iiasa/saved_vars/'
except:
	wddata='C:/Users/garyk/Documents/python_code/riskAssessmentFromPovertyEstimations/data/'
	wdfigs='C:/Users/garyk/Documents/python_code/riskAssessmentFromPovertyEstimations/figs/'
	wdvars='C:/Users/garyk/Documents/python_code/riskAssessmentFromPovertyEstimations/vars/'

MakePlots=False

latsubsaharan=np.load(wdvars+'latsubsaharan.npy')
lonsubsaharan=np.load(wdvars+'lonsubsaharan.npy')
africaMask1=np.load(wdvars+'africaMasksubsaharan.npy')

nyears=16

######################################
# Countries
######################################
	
f=open(wddata+'boundaries/countries_countryCodes.csv')
code=np.zeros(shape=(247),dtype=int)
countryNames=[]
i=-1
for line in f:
	i+=1
	tmp=line.split(',')
	code[i]=int(tmp[3])
	countryNames.append(tmp[0])

f=open(wddata+'boundaries/africanCountries.csv','r')
africanCountries=[]
for line in f:
	africanCountries.append(line[:-1])

#indexing the country codes
countryToIndex={}
indexToCountry={}
indexedcodes=np.zeros(shape=len(africanCountries))
for i in range(len(africanCountries)):
	j=np.where(africanCountries[i]==np.array(countryNames))[0][0]
	indexedcodes[i]=code[j]
	countryToIndex[africanCountries[i]]=indexedcodes[i]
	indexToCountry[indexedcodes[i]]=africanCountries[i]

countryToi={}
iToCountry={}
for i in range(len(indexedcodes)):
	index=indexedcodes[i]
	country=indexToCountry[index]
	countryToi[country]=i
	iToCountry[i]=country

######################################
# Sum each country
######################################
malCountry = np.load(wdvars+'malCountry_2000-2015.npy')

if MakePlots:
	year=1999
	for y in range(nyears):
		year+=1
		plt.clf()
		plt.imshow(malCountryGrid[y],cmap=cm.jet,vmax=20)
		plt.colorbar()
		plt.yticks([])
		plt.xticks([])
		plt.title(str(year)+': Average Malnutrition Prevalence')
		plt.savefig(wdfigs+'malCountryGrid_'+str(year)+'.png',dpi=700)

lats=np.zeros(shape=(44,16))
f = open(wddata+'boundaries/countryLats.csv')
k=-1
for line in f:
	line=line[:-1]
	tmp=line.split(',')
	k+=1
	if np.amax(tmp[3]==np.array(africanCountries))==0:
		continue
	i=countryToi[tmp[3]]
	lats[i,:]=float(tmp[1])
latsS=scale(lats)

######################################
# UNICEF numbers
######################################

MAMcaseload=np.zeros(shape=(44,nyears))

with open(wddata+'unicef_caseload/UNICEF_Wasting_2018_May/Trend-Table1.csv', 'rb') as csvfile:
	k=0
	reader = csv.reader(csvfile, delimiter=',', quotechar='"')
	for tmp in reader:
		k+=1

		country=tmp[2]
		if np.amax(country==np.array(africanCountries))==0:
			continue
		year=int(tmp[5])
		y=year-2000
		if year<2000 or year>2015:
			continue
		i=countryToi[country]
		MAMcaseload[i,y]=tmp[9]

MAMcaseload=np.ma.masked_array(MAMcaseload,MAMcaseload==0)
MAMmask=np.ma.getmask(MAMcaseload)
malCountryM=np.ma.masked_array(malCountry,MAMmask)
CorrM=corr(np.ma.compressed(MAMcaseload),np.ma.compressed(malCountryM))


SAMcaseload=np.zeros(shape=(44,nyears))

with open(wddata+'unicef_caseload/UNICEF_Severe_Wasting_2018_May/Trend-Table1.csv', 'rb') as csvfile:
	k=0
	reader = csv.reader(csvfile, delimiter=',', quotechar='"')
	for tmp in reader:
		k+=1

		country=tmp[2]
		if np.amax(country==np.array(africanCountries))==0:
			continue
		year=int(tmp[5])
		y=year-2000
		if year<2000 or year>2015:
			continue
		i=countryToi[country]
		SAMcaseload[i,y]=tmp[9]

SAMcaseload=np.ma.masked_array(SAMcaseload,SAMcaseload==0)
SAMmask=np.ma.getmask(SAMcaseload)
malCountryM=np.ma.masked_array(malCountry,SAMmask)
CorrS=corr(np.ma.compressed(SAMcaseload),np.ma.compressed(malCountryM))

#### MAM ####
x=np.ma.compressed(MAMcaseload)
ydata=np.ma.compressed(np.ma.masked_array(malCountry,MAMmask))
colors=np.ma.compressed(np.ma.masked_array(lats,MAMmask))

slope,b=np.polyfit(x,ydata,1)
yfit=slope*x+b
stdDev=np.std(ydata-x)

fig = plt.figure(figsize=(9, 5))
plt.clf()
plt.scatter(x,ydata,c=colors,marker='*',cmap=cm.jet)
plt.colorbar()
plt.plot(x,yfit,'-g')
plt.plot(x,yfit-stdDev,'--g',linewidth=0.5)
plt.plot(x,yfit+stdDev,'--g',linewidth=0.5)
plt.title('MAM vs Interpolated, Corr = '+str(round(CorrM,2))+', stdDev = '+str(round(stdDev,2)))
plt.xlabel('MAM %')
plt.ylabel('Interpolated %')
plt.savefig(wdfigs+'MAMvsInterpolated.pdf')
#############
#### SAM ####
x=np.ma.compressed(SAMcaseload)
ydata=np.ma.compressed(np.ma.masked_array(malCountry,SAMmask))
colors=np.ma.compressed(np.ma.masked_array(lats,SAMmask))

a,b,c=np.polyfit(x,ydata,2)
yfit=a*x**2+b*x+c
stdDev=np.std(ydata-x)
s=np.argsort(yfit)
yfit=yfit[s]
x2=x[s]

plt.clf()
plt.scatter(x,ydata,c=colors,marker='*',cmap=cm.jet)
plt.colorbar()
plt.plot(x2,yfit,'-g')
plt.plot(x2,yfit-stdDev,'--g',linewidth=0.5)
plt.plot(x2,yfit+stdDev,'--g',linewidth=0.5)
plt.title('SAM vs Interpolated, Corr = '+str(round(CorrS,2))+', stdDev = '+str(round(stdDev,2)))
plt.xlabel('SAM %')
plt.ylabel('Interpolated %')
plt.savefig(wdfigs+'SAMvsInterpolated.pdf')
#############




