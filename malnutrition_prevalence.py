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
import scipy
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
from libtiff import TIFF
import googlemaps

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

def createMidpointGridlatlon(grid,pixelsize):
	gridMid=np.zeros(shape=(grid.shape))
	for ilat in range(len(grid[:,0,0])):
		gridMid[ilat,:,0]=np.mean([grid[ilat,0,0],grid[ilat,0,0]+pixelsize])
	for ilon in range(len(grid[0,:,1])):
		gridMid[:,ilon,1]=np.mean([grid[:,ilon,1],grid[:,ilon,1]+pixelsize])
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
	closestDistM=np.swapaxes(closestDistM,0,1)
	iclosestDistM=np.swapaxes(iclosestDistM,0,1)
	return closestDistM,iclosestDistM

def nearestCity(cityLonLat,zeroLonLat):
	'''Give it:
	1. lon lat data set cities
	2. lon lat data set zeros on travel time
	Will return:
	1. Closest dist to city
	2. Closest city
	'''
	closestDist=10000*np.ones(shape=(len(zeroLonLat)))
	iclosestDist=np.ones(shape=(len(zeroLonLat)),dtype=int)
	R = 6373.0 #earth's radius
	for i in range(len(zeroLonLat)):
		distances=np.sqrt((zeroLonLat[i,0]-cityLonLat[:,0])**2+(zeroLonLat[i,1]-cityLonLat[:,1])**2)
		leastDist=np.amin(distances)
		iclosestDist[i]=np.where(distances==leastDist)[0][0]
		closestDist[i]=geopy.distance.distance([zeroLonLat[i,1],zeroLonLat[i,0]],[cityLonLat[iclosestDist[i],1],cityLonLat[iclosestDist[i],0]]).km
		print np.round(100*i/float(len(zeroLonLat)),2),'%' 

	return closestDist,iclosestDist

################################################
try:
	wddata='/Users/lilllianpetersen/iiasa/data/'
	wdfigs='/Users/lilllianpetersen/iiasa/figs/'
	wdvars='/Users/lilllianpetersen/iiasa/saved_vars/'
	tifWasting=TIFF.open(wddata+'malnutrition/IHME_AFRICA_CGF_2000_2015_WASTING_MEAN_2010_PREVALENCE_Y2018M02D28.TIF',mode='r')
except:
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
pixelsizem=abs(gt[-1])

# lat and lon
latm=np.ones(shape=(height1))
lonm=np.ones(shape=(width1))
for w in range(width1):
	lonm[w]=minx+w*pixelsizem
for h in range(height1):
	latm[h]=miny+h*pixelsizem
latm=latm[::-1] # reverse the order

mal=tifWasting.read_image()
mal[mal>1]=0.0308212 # The average of this weird pixel's neighbors

mal=mal[:1740]
latm=latm[:1740]
mal=mal[:,185:]
lonm=lonm[185:]

gridm=np.zeros(shape=(len(lonm),len(latm),2))
for x in range(len(lonm)):
	gridm[x,:,0]=lonm[x]
for y in range(len(latm)):
	gridm[:,y,1]=latm[y]

africaMask=np.zeros(shape=(mal.shape))
africaMask[mal<0]=1

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


######################################
# Travel time from 250k
######################################

ds=gdal.Open(wddata+'travel_time/TT_250K--SSA.tif')
width = ds.RasterXSize
height = ds.RasterYSize
gt = ds.GetGeoTransform()
minx = gt[0]
miny = gt[3] + width*gt[4] + height*gt[5] 
maxx = gt[0] + width*gt[1] + height*gt[2]
maxy = gt[3] 
pixelsize=abs(gt[-1])

latc=np.ones(shape=(height))
lonc=np.ones(shape=(width))
for w in range(width):
	lonc[w]=minx+w*pixelsize
for h in range(height):
	latc[h]=miny+h*pixelsize
latc=latc[::-1]
lattravel=latc
lontravel=lonc

travel=ds.ReadAsArray()
travel[travel<0]=-10

africaMask1=np.zeros(shape=(travel.shape),dtype=int)
africaMask1[travel<0]=1

######### 5km #########
latsubsaharan=latm[latm<np.amax(latc)+pixelsize]
latsubsaharan=latsubsaharan[latsubsaharan>np.amin(latc)]
travel=travel[latc>latm[-1]]
africaMask1=africaMask1[latc>latm[-1]]
latc=latc[latc>latm[-1]]
latsubsaharan=latsubsaharan[latsubsaharan>latm[-1]]

# travel
latl=np.radians(latc[::-1])+1.2
lonl=np.radians(lonc)+1.2
lut=RectSphereBivariateSpline(latl, lonl, travel)
newLats,newLons=np.meshgrid(np.radians(latsubsaharan[::-1])+1.2,np.radians(lonm)+1.2)
travel=lut.ev(newLats.ravel(),newLons.ravel()).reshape((len(lonm),len(latsubsaharan))).T

# AfricaMask
lut=RectSphereBivariateSpline(latl, lonl, africaMask1)
newLats,newLons=np.meshgrid(np.radians(latsubsaharan[::-1])+1.2,np.radians(lonm)+1.2)
africaMask1=lut.ev(newLats.ravel(),newLons.ravel()).reshape((len(lonm),len(latsubsaharan))).T
#######################
africaMask1=np.round(africaMask1,0)
lonsubsaharan=lonm

travel=np.ma.masked_array(travel,africaMask1)

plt.clf()
plt.imshow(travel,cmap=cm.gist_ncar_r,vmin=0,vmax=30)
plt.title('Travel Time')
plt.yticks([])
plt.xticks([])
plt.colorbar()
plt.savefig(wdfigs+'travel',dpi=700)


######################################
# Major Cities
######################################

#majorcities=np.zeros(shape=(210,2))
#majorCityNames=[]
#majorCityCountries=[]
#majorCityPop=np.zeros(shape=(210))
#f=open(wddata+'travel_time/african_major_cities.csv','r')
#i=-1
#for line in f:
#	i+=1
#	tmp=line.split(',')
#	majorCityNames.append(tmp[1])
#	majorCityCountries.append(tmp[2])
#	majorCityPop[i]=tmp[3]
#	if majorCityNames[i]=='Mbeya':
#	    majorcities[i]=-8.91667,33.4166
#	    continue
#	if majorCityNames[i]=='Kolwezi':
#	    majorcities[i]=-10.75000333,25.50000332
#	    continue
#	if majorCityNames[i]=='El Mahallah el Kubra':
#		majorcities[i]=30.9697,31.1681
#	else:
#		majorcities[i]=geolocator.geocode(str(majorCityNames[i]+', '+majorCityCountries[i])).latitude,geolocator.geocode(str(majorCityNames[i]+', '+majorCityCountries[i])).longitude
#	print i,majorCityNames[i],majorCityCountries[i],majorcities[i,0],majorcities[i,1]
    

cities=np.array(travel)
cities=np.ma.masked_array(cities,africaMask1)
cities[cities>.5]=1
cities[cities<1]=0
cities=np.ma.masked_array(cities,africaMask1)

cityrad=np.zeros(shape=(cities.shape))
citycenters=np.zeros(shape=(cities.shape))
for ilat in range(len(cities[:,0])-1,-1,-1):
	print np.round(100*ilat/float(len(cities[:,0])),2),'%'
	for ilon in range(len(cities[:,1])-1,-1,-1):
		if cities[ilat,ilon]==0 and cityrad[ilat,ilon]==0:
			citycenters[ilat,ilon]=1
			cityrad[ilat-10:ilat+11,ilon-10:ilon+11]=1

# for i in citycenters:
#     for j in i:
#         if 
#     exit()
    
gridsubsaharan=np.zeros(shape=(len(latsubsaharan),len(lonsubsaharan),2))
for x in range(len(latsubsaharan)):
	gridsubsaharan[x,:,0]=latsubsaharan[x]
	# lonz.append(lonc[x])
for y in range(len(lonsubsaharan)):
	gridsubsaharan[:,y,1]=lonsubsaharan[y]

midpointsc=createMidpointGridlatlon(gridsubsaharan,pixelsize) # lon lat

citycenters=np.ma.masked_array(citycenters,africaMask1)
centersLatLon=gridsubsaharan[citycenters==1]

closestcity,iclosestcity=nearestCity(majorcities,centersLatLon)

IDofmarketsheds=[]
citymatches=[]
matchcountries=[]
for cityindex in iclosestcity:
    IDofmarketsheds.append(majorCityNames[cityindex])
    citymatches.append(majorcities[cityindex])
    matchcountries.append(majorCityCountries[cityindex])

#f=open(wddata+'travel_time/citymatches_coord.csv','w')
#for i in range(len(centersLatLon)):
#    f.write(str(IDofmarketsheds[i]) +','+ str(matchcountries[i]) +','+ str(centersLatLon[i,0]) + "," + str(centersLatLon[i,1])+'\n')
#f.close()

revisedcities=np.zeros(shape=(91,2))
IDofmarketsheds=[]
matchcountries=[]
f=open(wddata+'travel_time/citymatches_coord.csv','r')
i=-1
for line in f:
    i+=1
    tmp=line.split(',')
    IDofmarketsheds.append(tmp[0])
    matchcountries.append(tmp[1])
    revisedcities[i,0]=tmp[2]
    revisedcities[i,1]=tmp[3]

closestcity2,iclosestcity2=nearestCity(majorcities,revisedcities)
# findNearest2(centersLatLon,midpointsc)
# closestcity,iclosestcity=nearestCity(majorcities,centersLatLon) OOF


######################################
# African Marketsheds
######################################

centersLonLat=np.zeros(shape=(revisedcities.shape))
centersLonLat[:,0]=revisedcities[:,1]
centersLonLat[:,1]=revisedcities[:,0]
marketSheds=np.zeros(shape=(travel.shape))
africaMaskLonLat=np.swapaxes(africaMask1,0,1)
midpointscLonLat=np.zeros(shape=(midpointsc.shape))
midpointscLonLat[:,:,0]=midpointsc[:,:,1]
midpointscLonLat[:,:,1]=midpointsc[:,:,0]
midpointscLonLat=np.swapaxes(midpointscLonLat,0,1)

plt.clf()
plt.plot(centersLonLat[:,0],centersLonLat[:,1],'*')
plt.title('African Marketsheds')
plt.savefig(wdfigs+'cities',dpi=700)

shedsDist,isheds=findNearest(centersLonLat, midpointscLonLat, africaMaskLonLat)


####### Sum Demand #######
malnum=np.array(malnumber)
malnum=np.ma.masked_array(malnum,africaMask1)

malnum=malnum[latm<np.amax(latsubsaharan)+pixelsizem+0.0001]
latms=latm[latm<np.amax(latsubsaharan)+pixelsizem+0.0001]
malnum=malnum[latms>np.amin(latsubsaharan)]
latms=latms[latms>np.amin(latsubsaharan)]

pop=pop[:,lonp<np.amax(lonm)+pixelsize]
lonp=lonp[lonp<np.amax(lonm)+pixelsize]
pop=pop[:,lonp>np.amin(lonm)]
lonp=lonp[lonp>np.amin(lonm)]

#### Number of grammes per child
tonnes=malnum*(365/75.)*(200)*(1/1000000.)*50

plt.clf()
plt.imshow(tonnes,cmap=cm.jet,vmin=0,vmax=10)
plt.title('Tonnes of SC+ Required per Year')
plt.yticks([])
plt.xticks([])
plt.colorbar()
plt.savefig(wdfigs+'tonnesSCperyear',dpi=700)
################################

shedsDemand=np.zeros(shape=(np.amax(isheds)+1))
ishedsDemand=np.zeros(shape=(isheds.shape))
ishedsDemandScaled=np.zeros(shape=(isheds.shape))
for i in range(np.amax(isheds)+1):
	shedsDemand[i]=np.sum(tonnes[isheds==i])
	ishedsDemand[isheds==i]=shedsDemand[i]
	ishedsDemandScaled[isheds==i]=shedsDemand[i]/(float(np.where(isheds==i)[1].shape[0]))

ishedsDemand=np.ma.masked_array(ishedsDemand,africaMask1)
ishedsDemandScaled=np.ma.masked_array(ishedsDemandScaled,africaMask1)

plt.clf()
plt.imshow(ishedsDemand,cmap=cm.jet,vmin=0,vmax=20000)
plt.title('Tonnes SC+ Required by Marketshed')
plt.yticks([])
plt.xticks([])
plt.colorbar()
plt.savefig(wdfigs+'ishedsDemand',dpi=700)

plt.clf()
plt.imshow(malnum,cmap=cm.jet,vmin=0,vmax=100)
plt.title('Number of Malnurished')
plt.yticks([])
plt.xticks([])
plt.colorbar()
plt.savefig(wdfigs+'malnumsubsaharan',dpi=700)


plt.clf()
plt.imshow(isheds,cmap=cm.nipy_spectral)
plt.title('Market Sheds')
plt.colorbar()
plt.savefig(wdfigs+'marketsheds',dpi=700)

plt.clf()
plt.imshow(shedsDist,cmap=cm.nipy_spectral)
plt.title('Dist to Market Sheds')
plt.colorbar()
plt.savefig(wdfigs+'marketshedsDist',dpi=700)


plt.clf()
plt.imshow(citycenters,cmap=cm.nipy_spectral_r)
plt.yticks([])
plt.xticks([])
plt.title('travel time to 250k')
plt.colorbar()
plt.savefig(wdfigs +'citiescenters',dpi=900)

cities=np.array()

plt.clf()
plt.imshow(travel,cmap=cm.nipy_spectral,vmin=1,vmax=2)
plt.yticks([])
plt.xticks([])
plt.title('travel time to 250k')
plt.colorbar()
plt.savefig(wdfigs +'traveltime250k',dpi=700)


######################################
# mapping between
######################################
exit()
gmaps = googlemaps.Client(key='AIzaSyAv4HITl2PsxqID8CX8xbOa8qMv6CU03hA')
distanceArray=np.zeros(shape=(91,91))
distanceDictionary={}
counter=0
listofcities=IDofmarketsheds
for i in range(len(listofcities)):
    print(listofcities[i])
    distanceDictionary[listofcities[i]]=[]
    for j in range(len(listofcities)):
        if listofcities[i]==listofcities[j]:
            distanceDictionary[listofcities[i]].append(0)
            distanceArray[i,j]=0
        else:
            gmapreturn=(gmaps.distance_matrix(listofcities[i]+" "+matchcountries[i],listofcities[j]+" "+matchcountries[j])['rows'][0]['elements'][0])
            if(gmapreturn=={u'status': u'ZERO_RESULTS'} or gmapreturn=={u'status': u'NOT_FOUND'}):
                distanceDictionary[listofcities[i]].append(999999999)
                distanceArray[i,j]=999999999
                print('wrong!')
            else:
                distanceDictionary[listofcities[i]].append(gmapreturn['distance']['value'])
                distanceArray[i,j]=gmapreturn['distance']['value']

for key in distanceDictionary:
    for i in range(len(distanceDictionary[key])):
        distanceDictionary[key][i]=(distanceDictionary[key][i]/1000.0)

distanceArray=distanceArray/1000.0

##9cents per metric tonne km
transportcostArray=distanceArray*0.09

networkDictionary={}
### array to dict
for i in range(len(listofcities)):
    networkDictionary[listofcities[i]]=[]
    for j in range(len(listofcities)):
        networkDictionary[listofcities[i]].append(transportcostArray[i,j])

with open(wddata + 'travel_time/traveldictionary.json', 'w') as outfile:
      json.dump(distanceDictionary, outfile)

with open(wddata + 'travel_time/traveldictionary2.json', 'w') as outfile:
      json.dump(networkDictionary, outfile)

df = pd.DataFrame(distanceArray)
df.to_csv(wddata+"travelarray.csv", header=None, index=None)

with open(wddata + 'travel_time/traveldictionary.json') as f:
      distanceDictionary=json.load(f)

df = pd.DataFrame(distanceArray)
df.to_csv(wddata+"travelarray.csv", header=None, index=None)

distanceArray = np.genfromtxt(wddata+'travelarray.csv', delimiter=',')

###convert to cost of transport per tonne
# distanceDictionary.update((x, y/1000*) for x, y in distanceDictionary.items())
subsaharancountry=[]
f=open(wddata+'population/subsaharan.csv','r')
i=-1
for line in f:
    i+=1
    tmp=line.split(',')
    subsaharancountry.append(tmp[0][:-1])

countryW=[]
wastingnational=np.zeros(shape=193)
SAMnational=np.zeros(shape=193)
stuntingnational=np.zeros(shape=193)
f=open(wddata+'population/nationalwastingnumbers.csv','r')
i=-1
for line in f:
    i+=1
    tmp=line.split(',')
    countryW.append(tmp[0])
    try:
        wastingnational[i]=float(tmp[1])
    except:
        wastingnational[i]=0
    try:
        SAMnational[i]=float(tmp[2])
    except:
        SAMnational[i]=0
    try:
        stuntingnational[i]=float(tmp[3])
    except:
        stuntingnational[i]=0
  
indexedwasting=np.zeros(shape=50)
indexedSAM=np.zeros(shape=50)
indexedstunting=np.zeros(shape=50)
indexedMAM=np.zeros(shape=50)
for j in range(len(subsaharancountry)):
    for i in range(len(countryW)):
        if(subsaharancountry[j] == countryW[i]):
            indexedwasting[j]=wastingnational[i]*1000.0
            indexedSAM[j]=SAMnational[i]*1000.0
            indexedstunting[j]=stuntingnational[i]*1000.0
            indexedMAM[j]=indexedwasting[j]-indexedSAM[j]

f=open(wddata+'population/nationalcasenumbers.csv','w')
for i in range(len(i)):
    f.write(str(subsaharancountry[i]) +','+ str(indexedwasting[i]) +','+ str(indexedSAM[i]) + "," + str(indexedMAM[i]) + "," + str(indexedstunting[i])+'\n')
f.close()