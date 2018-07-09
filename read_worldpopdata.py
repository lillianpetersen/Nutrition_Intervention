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
from math import sin, cos, sqrt, atan2, radians
from geopy.geocoders import Nominatim
geolocator = Nominatim()
import geopy.distance

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

from itertools import compress
class MaskableList(list):
    def __getitem__(self, index):
        try: return super(MaskableList, self).__getitem__(index)
        except TypeError: return MaskableList(compress(self, index))

def checkIfInPolygon(lon, lat,polygon):
	point = Point(lon, lat)
	return polygon.contains(point)

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
		if yClosestL<midPoint[1]:
			yClosestGrid=np.where(grid[0,:,1]==yClosestL)[0][0]
		else:
			yClosestGrid=np.where(grid[0,:,1]==yClosestL)[0][0]+1
		if xClosestL<midPoint[0]:
			xClosestGrid=np.where(grid[:,0,0]==xClosestL)[0][0]
		else:
			xClosestGrid=np.where(grid[:,0,0]==xClosestL)[0][0]-1
		midpointHits[i]=xClosestGrid,yClosestGrid

		### Convert radians to meters ###
		R = 6373.0 # Earth's radius
		lon1,lon2=radians(x1),radians(x2)
		lat1,lat2=radians(y1),radians(y2)
		dlon = lon2 - lon1
		dlat = lat2 - lat1
		a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
		c = 2 * atan2(sqrt(a), sqrt(1 - a))

		dist[i] = R * c
		
	return midpointHits,dist

def createMidpointGrid(grid,pixelsize):
	gridMid=np.zeros(shape=(grid.shape))
	for ilon in range(len(grid[:,0,0])):
		gridMid[ilon,:,0]=np.mean([grid[ilon,0,0],grid[ilon,0,0]+pixelsize])
	for ilat in range(len(grid[0,:,1])):
		gridMid[:,ilat,1]=np.mean([grid[:,ilat,1],grid[:,ilat,1]+pixelsize])
	return gridMid

def findNearestCities(cityLonLat,gridMid,imageMask1):
closestDist=10000*np.ones(shape=(gridMid[:,:,0].shape))
closestDistDeg=100*np.ones(shape=(gridMid[:,:,0].shape))
iclosestDist=np.zeros(shape=(gridMid[:,:,0].shape))
imageMask4=imageMask1
R = 6373.0 #earth's radius
for ilon in range(len(gridMid[:,0,0])):
	for ilat in range(len(gridMid[0,:,1])):
		if imageMask1[ilon,ilat]==True:
			continue
		lon,lat=gridMid[ilon,0,0],gridMid[0,ilat,1]
		for city in range(len(cityLonLat[:,0])):
			clon,clat=cityLonLat[city,0],cityLonLat[city,1]
			dist=sqrt((clon-lon)**2+(clat-lat)**2)
			if dist<closestDistDeg[ilon,ilat]:
				closestDistDeg[ilon,ilat]=dist
				iclosestDist[ilon,ilat]=city
				clonBest,clatBest=clon,clat
		
		coord1=[lat,lon]
		coord2=[clatBest,clonBest]
		closestDist[ilon,ilat]=geopy.distance.vincenty(coord1,coord2).km
	print np.round(100*ilon/float(len(gridMid[:,0,0])),2),'%'
closestDistM=np.ma.masked_array(closestDist,imageMask1)
closestDistDegM=np.ma.masked_array(closestDistDeg,imageMask1)
iclosestDistM=np.ma.masked_array(iclosestDist,imageMask1)
closestDistM=np.swapaxes(closestDistM,0,1)
closestDistDegM=np.swapaxes(closestDistDegM,0,1)
iclosestDistM=np.swapaxes(iclosestDistM,0,1)
		return closestDistM,iclosestDistM


	

###############################################

wddata='/Users/lilllianpetersen/iiasa/data/'
wdfigs='/Users/lilllianpetersen/iiasa/figs/'
wdvars='/Users/lilllianpetersen/iiasa/saved_vars/'


countriesT=['Malawi','Nigeria','Uganda','Tanzania']
countries=[]
for c in countriesT:
	countries.append(c.lower())
countryAbb=['mwi','nga','uga','tza']
years=['11','10','10','10']

#i=-1
#for country in countries:
#	i+=1
i=1 #country=Nigeria
country='Nigeria'


tif125 = TIFF.open(wddata+'poverty/'+country+'/'+countryAbb[i]+years[i]+'povcons125.tif',mode='r')
ds=gdal.Open(wddata+'poverty/'+country+'/'+countryAbb[i]+years[i]+'povcons125.tif')
width = ds.RasterXSize
height = ds.RasterYSize
gt = ds.GetGeoTransform()
minx = gt[0]
miny = gt[3] + width*gt[4] + height*gt[5] 
maxx = gt[0] + width*gt[1] + height*gt[2]
maxy = gt[3] 
pixelsize=abs(gt[-1])

lat=np.ones(shape=(height))
lon=np.ones(shape=(width))
for w in range(width):
	lon[w]=minx+w*pixelsize
for h in range(height):
	lat[h]=miny+h*pixelsize
lat=lat[::-1] # reverse the order

# read your shapefile
plt.clf()

tif125 = TIFF.open(wddata+'poverty/'+country+'/'+countryAbb[i]+years[i]+'povcons125.tif',mode='r')
tif200 = TIFF.open(wddata+'poverty/'+country+'/'+countryAbb[i]+years[i]+'povcons200.tif',mode='r')
tifuncert200 = TIFF.open(wddata+'poverty/'+country+'/'+countryAbb[i]+years[i]+'povcons200uncert.tif',mode='r')
tifuncert125 = TIFF.open(wddata+'poverty/'+country+'/'+countryAbb[i]+years[i]+'povcons125uncert.tif',mode='r')

pov125 = tif125.read_image()*100
pov200 = tif200.read_image()*100
uncert125 = tifuncert125.read_image()*100
uncert200 = tifuncert200.read_image()*100

imageMask=pov125<0

pov125=np.ma.masked_array(pov125,imageMask)
pov125=pov125[~imageMask.all(axis=1)]
pov125=pov125[:,~imageMask.all(axis=0)]

lonM=np.ma.masked_array(lon,imageMask.all(axis=0))
lonM=np.ma.compressed(lonM)
latM=np.ma.masked_array(lat,imageMask.all(axis=1))
latM=np.ma.compressed(latM)

pov200=np.ma.masked_array(pov200,imageMask)
pov200=pov200[~imageMask.all(axis=1)]
pov200=pov200[:,~imageMask.all(axis=0)]

uncert125=np.ma.masked_array(uncert125,imageMask)
uncert125=uncert125[~imageMask.all(axis=1)]
uncert125=uncert125[:,~imageMask.all(axis=0)]

uncert200=np.ma.masked_array(uncert200,imageMask)
uncert200=uncert200[~imageMask.all(axis=1)]
uncert200=uncert200[:,~imageMask.all(axis=0)]

plt.clf()
plt.imshow(pov125,cmap=cm.hot_r,vmin=0,vmax=100)
plt.yticks([])
plt.xticks([])
plt.title('Percent living under 1.25/day '+countriesT[i]+' 20'+years[i])
plt.colorbar()
plt.savefig(wdfigs+country+'_poverty_125',dpi=700)

#plt.clf()
#plt.imshow(pov200,cmap=cm.hot_r,vmin=0,vmax=100)
#plt.yticks([])
#plt.xticks([])
#plt.title('Percent living under 2USD/day '+countriesT[i]+' 20'+years[i])
#plt.colorbar()
#plt.savefig(wdfigs+country+'_poverty_200',dpi=700)
#
#plt.clf()
#plt.imshow(uncert200,cmap=cm.winter_r,vmin=0,vmax=100)
#plt.yticks([])
#plt.xticks([])
#plt.title(countriesT[i]+' 95 Percent Confidence Interval for Pov2USD')
#plt.colorbar()
#plt.savefig(wdfigs+country+'_uncertainty200',dpi=700)
#
#plt.clf()
#plt.imshow(uncert125,cmap=cm.winter_r,vmin=0,vmax=100)
#plt.yticks([])
#plt.xticks([])
#plt.title(countriesT[i]+' 95 Percent Confidence Interval for Pov1.25USD')
#plt.colorbar()
#plt.savefig(wdfigs+country+'_uncertainty125',dpi=700)

######################################################
# Roads
######################################################

country='nigeria'
rRoads = shapefile.Reader(wddata+'openstreetmap/'+country+'/openstreetmap/nga_trs_roads_osm.shp')
#rRoads = shapefile.Reader(wddata+'nigeria_roads/part/nga_trs_roads_osm_wfpschema_sep2016/nga_trs_roads_osm_wfpschema_sep2016.shp')

records=rRoads.shapeRecords()
records=MaskableList(records)
roadType=[rec.record[6] for rec in records]

maskP=np.array([rType=='primary' for rType in roadType])
maskS=np.array([rType=='secondary' for rType in roadType])
maskT=np.array([rType=='tertiary' for rType in roadType])

records1=records[maskP]
records2=records[maskS]
records3=records[maskT]


#plt.clf()
#map = Basemap(llcrnrlon=lonM[0],llcrnrlat=np.amin(latM),urcrnrlon=lonM[-1],urcrnrlat=np.amax(latM), projection='lcc',lon_0=(lonM[-1]+lonM[0])/2,lat_0=(latM[-1]+latM[0])/2,resolution='i')
#map.drawcoastlines(linewidth=1.5)
#map.drawcountries()
#print 'countries'
#
#map.readshapefile(wddata+'openstreetmap/nigeria/openstreetmap/nga_trs_roads_osm', name='roads', drawbounds=True)
#
#ax = plt.gca()
#plt.title('Nigerian Roads from Open Street Map')
#
#plt.savefig(wdfigs+'openstreetmap_nigeria',dpi=700)

grid = np.zeros(shape=(len(lonM),len(latM),2))
for x in range(len(lonM)):
	grid[x,:,0]=lonM[x]
for y in range(len(latM)):
	grid[:,y,1]=latM[y]

try:
	roadDensity=np.load(wdvars+'nigeria_raodDensity.npy')
except:
	roadDensity=np.zeros(shape=(len(lonM),len(latM),3))
	
	for i in range(len(records1)):
		print i,np.round(100*i/float(len(records1)),1)
		shapePoints=records1[i].shape.points
		midpointHits,dist=findGridAndDistance(shapePoints,grid)
		midpointHits=np.array(midpointHits,dtype=int)
		for j in range(len(midpointHits)):
			roadDensity[midpointHits[j,0],midpointHits[j,1],0]+=dist[j]
	
	for i in range(len(records2)):
		print i,np.round(100*i/float(len(records2)),1)
		shapePoints=records2[i].shape.points
		midpointHits,dist=findGridAndDistance(shapePoints,grid)
		midpointHits=np.array(midpointHits,dtype=int)
		for j in range(len(midpointHits)):
			roadDensity[midpointHits[j,0],midpointHits[j,1],1]+=dist[j]
	
	for i in range(len(records3)):
		print i,np.round(100*i/float(len(records3)),2)
		shapePoints=records3[i].shape.points
		midpointHits,dist=findGridAndDistance(shapePoints,grid)
		midpointHits=np.array(midpointHits,dtype=int)
		for j in range(len(midpointHits)):
			roadDensity[midpointHits[j,0],midpointHits[j,1],2]+=dist[j]
	
	np.save(wdvars+'nigeria_raodDensity',roadDensity)

imageMask2=imageMask
imageMask2=imageMask2[~imageMask2.all(axis=1)]
imageMask2=imageMask2[:,~imageMask2.all(axis=0)]
imageMask1=np.swapaxes(imageMask2,0,1)
imageMask3=np.zeros(shape=(imageMask2.shape[0],imageMask2.shape[1],3),dtype=bool)
for i in range(3):
	imageMask3[:,:,i]=imageMask2

roadDensityM=np.swapaxes(roadDensity[:,:,:],0,1)
roadDensityM=np.ma.masked_array(roadDensityM,imageMask3)

#plt.clf()
##plt.imshow(roadDensity,cmap=cm.gist_heat_r)
#plt.imshow(roadDensityM[:,:,0],cmap=cm.nipy_spectral_r)
#plt.colorbar()
#plt.title('km primary roads per square km in Nigeria')
#plt.savefig(wdfigs+'nigeria_primary_road_density',dpi=800)
#
#plt.clf()
#plt.figure(1,figsize=(3,4))
##plt.imshow(roadDensity,cmap=cm.gist_heat_r)
#plt.imshow(roadDensityM[:,:,1],cmap=cm.nipy_spectral_r)
#plt.colorbar()
#plt.title('Nigeria: km secondary roads/km^2')
#plt.savefig(wdfigs+'nigeria_secondary_road_density',dpi=800)
#
#plt.clf()
#plt.figure(1,figsize=(3,4))
#plt.imshow(roadDensityM[:,:,2],cmap=cm.nipy_spectral_r)
#plt.colorbar()
#plt.title('Nigeria: km tertiary roads/km^2')
#plt.savefig(wdfigs+'nigeria_tertiary_road_density',dpi=800)
#
#roads0=roadDensity>0
#roads0=1*roads0
#
#plt.clf()
#plt.figure(1,figsize=(3,4))
#plt.imshow(roads0[:,:,0],cmap=cm.gist_heat_r)
#plt.colorbar()
#plt.title('Primary Roads Nigeria')
#plt.savefig(wdfigs+'nigeria_primary_road_mask',dpi=600)
#
#plt.clf()
#plt.figure(1,figsize=(3,4))
#plt.imshow(roads0[:,:,1],cmap=cm.gist_heat_r)
#plt.colorbar()
#plt.title('Secondary Roads Nigeria')
#plt.savefig(wdfigs+'nigeria_secondary_road_mask',dpi=600)
#
#plt.clf()
##plt.imshow(roadDensity,cmap=cm.gist_heat_r)
#plt.imshow(roadDensityM[:,::-1],cmap=cm.nipy_spectral_r)
#plt.colorbar()
#plt.xticks([])
#plt.yticks([])
#plt.title('Nigerian Roads')
#plt.savefig(wdfigs+'nigeria_all_roads',dpi=800)


############### Claculate Distance to Roads ###############
try:
	distToRoads=np.load(wdvars+'nigeria_distToRoads.npy')
except:
	distToRoads=1-roads0 # road = 0
	distToRoads[distToRoads==1]=9999.
	
	for iroad in range(1,3):
		for i in range(300):
			print np.round(100*i/float(300),2),'%'
			for ilon in range(len(lonM)):
				if ilon==0:
					continue
				for ilat in range(len(latM)):
					if ilat==0:
						continue
					if distToRoads[ilon,ilat,iroad]==9999. and np.amin(distToRoads[ilon-1:ilon+2, ilat-1:ilat+2, iroad])==i:
						distToRoads[ilon,ilat,iroad]=i+1
		
	np.save(wdvars+'nigeria_distToRoads',distToRoads)
distToRoadsMask=distToRoads>=9999.
distToRoads=np.ma.masked_array(distToRoads,distToRoadsMask)
distToRoadstmp=np.swapaxes(distToRoads[:,:,:],0,1)
distToRoadsM=np.ma.masked_array(distToRoadstmp,imageMask3)

#print 'here'
#plt.clf()
#x,y = np.ma.compressed(grid[:,:,0]),np.ma.compressed(grid[:,:,1])
#plt.pcolor(x,y,distToRoadsM[:,:,0],cmap=cm.gist_heat_r)
##plt.imshow(distToRoadsM[:,:,0],cmap=cm.gist_heat_r)
#plt.colorbar()
#plt.xticks([])
#plt.yticks([])
#plt.scatter(cityLonLat[:,0], cityLonLat[:,1],s=popScaled,c=cityPop,cmap=cm.jet,alpha=0.7,vmin=0,vmax=800000)
#plt.title('Kilometers to Primary Roads Nigeria')
#plt.savefig(wdfigs+'nigeria_dist_to_primary_road',dpi=800)
#exit()

exit()
pov125M=np.ma.masked_array(pov125,np.swapaxes(distToRoadsMask[:,:,0],0,1))
var=['Primary','Secondary','Tertiary']
Corr=np.zeros(shape=(6))
slope=np.zeros(shape=(6))
bInt=np.zeros(shape=(6))
xMulti=np.zeros(shape=(1079449,3))
for i in range(3):
	x=np.ma.compressed(distToRoadsM[:,:,i])
	xMulti[:,i]=x
	ydata=np.ma.compressed(pov125M)
	yMulti=ydata
	Corr[i]=corr(x,ydata)
	slope[i],bInt[i]=np.polyfit(x,ydata,1)
	yfit=slope[i]*x+bInt[i]
	plt.clf()
	plt.figure(21,figsize=(8,6))
	plt.plot(x,ydata,'*',markersize=1)
	plt.plot(x,yfit,'g-')
	plt.xlabel('Dist to '+var[i]+' Roads')
	plt.ylabel('% Below 1.25$/day')
	plt.title('Distance to '+var[i]+' Roads and Poverty %, Corr = '+str(np.round(Corr[i],2)))
	plt.savefig(wdfigs+'roaddist_'+var[i]+'_pov_corr',dpi=700)

roadDensityM2=np.array(roadDensityM2)
roadDensityM2mask=roadDensityM==0
roadDensityM2=np.ma.masked_array(roadDensityM,roadDensityM2mask)
for i in range(3):
	pov125M=np.ma.masked_array(pov125,roadDensityM2mask[:,:,i])
	x=np.ma.compressed(roadDensityM2[:,:,i])
	print x.shape
	ydata=np.ma.compressed(pov125M)
	Corr[i+3]=corr(x,ydata)
	slope[i+3],bInt[i+3]=np.polyfit(x,ydata,1)
	yfit=slope[i+3]*x+bInt[i+3]
	plt.clf()
	plt.figure(21,figsize=(8,6))
	plt.plot(x,ydata,'*',markersize=1)
	plt.plot(x,yfit,'g-')
	plt.xlabel('Density of '+var[i]+' Roads')
	plt.ylabel('% Below 1.25$/day')
	plt.title('Density of '+var[i]+' Roads and Poverty %, Corr = '+str(np.round(Corr[i+3],2)))
	plt.savefig(wdfigs+'roaddensity_'+var[i]+'_pov_corr',dpi=700)

############## Multivariate ##############
#xMulti=np.zeros(shape=(k))
#cropYield1=np.ma.masked_array(cropYield,anomMask[:,:,3])
#ydata=np.ma.compressed(cropYield1)
#
#iMulti=-1
#for m in range(5):
#    iMulti+=1
#    x=np.ma.compressed(ndviAnom[:,:,m])
#    xMulti[:,iMulti]=x
#
#for m in range(5):
#    iMulti+=1
#    x=np.ma.compressed(eviAnom[:,:,m])
#    xMulti[:,iMulti]=x
#
#for m in range(5):
#    iMulti+=1
#    x=np.ma.compressed(ndwiAnom[:,:,m])
#    xMulti[:,iMulti]=x
from sklearn import linear_model

clf=linear_model.LinearRegression()
clf.fit(xMulti,yMulti)
exit()

np.save(wdvars+'Illinois/xMulti',xMulti)
np.save(wdvars+'Illinois/ydataMulti',ydata)


############### Claculate Distance to Cities ###############

try:
	cityPop=np.load('NigeriaCityPop.npy')
	cityLonLat=np.load('NigeriaCityLonLat.npy')
except:
	f = open(wddata+'openstreetmap/nigeria/nigeria_cities_pop.csv','r')
	
	cityPop=np.zeros(shape=(400))
	cityLonLat=np.zeros(shape=(400,2))
	city=[]
	i=-1
	for line in f:
		i+=1
		tmp=line.split(',')
		citytmp=tmp[0]
		city.append(citytmp)
		location = geolocator.geocode(citytmp+', Nigeria')
		cityLonLat[i]=location.longitude,location.latitude
		cityPop[i]=float(tmp[1])
	np.save(wdvars,'NigeriaCityPop',cityPop)
	np.save(wdvars,'NigeriaCityLonLat',cityLonLat)

try:
	closestDistM=np.load(wdvars+'closestDist_to_nigerian_cities.npy')
	iclosestDistM=np.load(wdvars+'iclosestDist_to_nigerian_cities.npy')
except:
	closestDistM,iclosestDistM=findNearestCities(cityLonLat,gridMid,imageMask1)
	np.save(wdvars+'closestDist_to_nigerian_cities',np.array(closestDistM))
	np.save(wdvars+'iclosestDist_to_nigerian_cities',np.array(iclosestDistM))

plt.clf()
plt.imshow(closestDistM,cmap=cm.viridis_r)
plt.colorbar()
plt.title('Distance to Nearest City')
plt.savefig(wdfigs+'closestDist',dpi=700)

dmin=0
dmax=np.amax(cityPop)
dif=dmax-dmin
popScaled=((cityPop-dmin)/dif)
popScaled=1-popScaled

popClosestDist=np.zeros(shape=(iclosestDistM.shape))
popClosestDistScaled=np.zeros(shape=(iclosestDistM.shape))
for ilat in range(len(iclosestDistM[:,0])):
	for ilon in range(len(iclosestDistM[0,:])):
		if imageMask2[ilat,ilon]==False:
			popClosestDistScaled[ilat,ilon]=popScaled[int(iclosestDistM[ilat,ilon])]
			popClosestDist[ilat,ilon]=cityPop[int(iclosestDistM[ilat,ilon])]
popClosestDistM=np.ma.masked_array(popClosestDist,imageMask2)

dmin=0
dmax=np.amax(closestDistM)
dif=dmax-dmin
max1=1
min1=.3
closestDistScaled=((closestDistM-dmin)/dif)*(max1-min1)+min1
distToCityIndex=closestDistScaled*popClosestDistScaled

plt.clf()
plt.imshow(distToCityIndex,cmap=cm.viridis_r,vmin=0,vmax=1)
plt.colorbar()
plt.title('Rural Index by Cities')
plt.savefig(wdfigs+'distToCityIndex',dpi=700)


plt.clf()
plt.imshow(iclosestDistM)
plt.colorbar()
plt.title('Nearest City')
plt.savefig(wdfigs+'iclosestDist',dpi=700)

plt.clf()
plt.imshow(popClosestDistM)
plt.colorbar()
plt.title('Nearest City Population')
plt.savefig(wdfigs+'popClosestDist',dpi=700)

#plt.clf()
#map = Basemap(llcrnrlon=lonM[0],llcrnrlat=np.amin(latM),urcrnrlon=lonM[-1],urcrnrlat=np.amax(latM), projection='lcc',lon_0=(lonM[-1]+lonM[0])/2,lat_0=(latM[-1]+latM[0])/2,resolution='i')
#map.drawcoastlines(linewidth=1.5)
#map.drawcountries()
#map.drawlsmask(land_color='w', ocean_color='lightblue',zorder=0)
#x,y = map(np.ma.compressed(grid[:,:,0]),np.ma.compressed(grid[:,:,1]))
#map.pcolormesh(x,y,roadDensity[:,:,0], latlon=False, cmap=cm.nipy_spectral_r,zorder=1)
#x,y = map(cityLonLat[:,0], cityLonLat[:,1])
#map.scatter(x,y,s=popScaled,c=cityPop,cmap=cm.jet,alpha=0.7,vmin=0,vmax=800000,zorder=2)
#plt.title('Nigerian Cities by Population')
#plt.colorbar()
#plt.savefig(wdfigs+'nigeria_cities',dpi=700)


exit()

######################################################################
# Age Structures
######################################################################
for y in range(0,66,5):
	for g in ('M','F'):
		if y<10:
			print 'A0'+str(y)+'0'+str(y+4)+g
			# load tiff
			vars()['tifA0'+str(y)+'0'+str(y+4)+g]=TIFF.open(wddata+'africa_1km_age_structures/AFR_PPP_A0'+str(y)+'0'+str(y+4)+'_'+g+'_2015_adj_v5.tif',mode='r')

			# convert tiff to array
			vars()['A0'+str(y)+'0'+str(y+4)+g]=vars()['tifA0'+str(y)+'0'+str(y+4)+g].read_image() 
			# Mask
			imageMask=vars()['A0'+str(y)+'0'+str(y+4)+g]<0
			vars()['A0'+str(y)+'0'+str(y+4)+g]=np.ma.masked_array(vars()['A0'+str(y)+'0'+str(y+4)+g],imageMask)
			vars()['A0'+str(y)+'0'+str(y+4)+g]=vars()['A0'+str(y)+'0'+str(y+4)+g][~imageMask.all(axis=1)]
			vars()['A0'+str(y)+'0'+str(y+4)+g]=vars()['A0'+str(y)+'0'+str(y+4)+g][:,~imageMask.all(axis=0)]
			
		elif y==65:
			print 'A'+str(y)+'PL'+g
			# load tiff
			vars()['tifA'+str(y)+'PL'+g]=TIFF.open(wddata+'africa_1km_age_structures/AFR_PPP_A'+str(y)+'PL'+'_'+g+'_2015_adj_v5.tif',mode='r')

			# convert tiff to array
			vars()['A'+str(y)+'PL'+g]=vars()['tifA'+str(y)+'PL'+g].read_image()

			imageMask=vars()['A'+str(y)+'PL'+g]<0
			vars()['A'+str(y)+'PL'+g]=np.ma.masked_array(vars()['A'+str(y)+'PL'+g],imageMask)
			vars()['A'+str(y)+'PL'+g]=vars()['A'+str(y)+'PL'+g][~imageMask.all(axis=1)]
			vars()['A'+str(y)+'PL'+g]=vars()['A'+str(y)+'PL'+g][:,~imageMask.all(axis=0)]
		else:
			print 'A'+str(y)+str(y+4)+g
			# load tiff
			vars()['tifA'+str(y)+str(y+4)+g]=TIFF.open(wddata+'africa_1km_age_structures/AFR_PPP_A'+str(y)+str(y+4)+'_'+g+'_2015_adj_v5.tif',mode='r')

			# convert tiff to array
			vars()['A'+str(y)+str(y+4)+g]=vars()['tifA'+str(y)+str(y+4)+g].read_image()

			imageMask=vars()['A'+str(y)+str(y+4)+g]<0
			vars()['A'+str(y)+str(y+4)+g]=np.ma.masked_array(vars()['A'+str(y)+str(y+4)+g],imageMask)
			vars()['A'+str(y)+str(y+4)+g]=vars()['A'+str(y)+str(y+4)+g][~imageMask.all(axis=1)]
			vars()['A'+str(y)+str(y+4)+g]=vars()['A'+str(y)+str(y+4)+g][:,~imageMask.all(axis=0)]

ageStructures=np.zeros(shape=(8,2))
gender=['M','F']
xc,yc=3738,1064
i=-1
for y in range(0,36,5):
	i+=1
	for g in range(2):
		if y<10:
			ageStructures[i,g]=vars()['A0'+str(y)+'0'+str(y+4)+gender[g]][xc,yc]
		else:
			ageStructures[i,g]=vars()['A'+str(y)+str(y+4)+gender[g]][xc,yc]

X = np.arange(0,36,5)


plt.clf()
plt.barh(X, ageStructures[:,1],height=4, color = 'm',alpha=.7,linewidth=2,edgecolor='k')
plt.barh(X, -ageStructures[:,0],height=4, color = 'b',alpha=.7,linewidth=2,edgecolor='k')
plt.axvline(x=0,color='k',linewidth=4)
plt.title('Age Pyramid in a 1km Grid Cell')
plt.ylabel('Ages')
plt.xlabel('Number')
plt.savefig(wdfigs+'sample_age_structure',dpi=700)

exit()


tifA0004F=TIFF.open(wddata+'africa_1km_age_structures/AFR_PPP_A0004_F_2015_adj_v5.tif',mode='r')
A0004F=tifA0004F.read_image()
tifA0004M=TIFF.open(wddata+'africa_1km_age_structures/AFR_PPP_A0004_M_2015_adj_v5.tif',mode='r')
A0004M=tifA0004F.read_image()
tifWOCBA=TIFF.open(wddata+'africa_1km_age_structures/AFR_PPP_WOCBA_2015_adj_v5.tif',mode='r')
WOCBA=tifWOCBA.read_image()

imageMask=A0004F<0

A0004F=np.ma.masked_array(A0004F,imageMask)
A0004M=np.ma.masked_array(A0004M,imageMask)
WOCBA=np.ma.masked_array(WOCBA,imageMask)

A0004=A0004F+A0004M

plt.clf()
plt.imshow(A0004,cmap=cm.jet,vmin=0,vmax=40)
plt.yticks([])
plt.xticks([])
plt.title('Number of Children Age 0-4 2015')
plt.colorbar()
plt.savefig(wdfigs+'A0004',dpi=700)

plt.clf()
plt.imshow(WOCBA,cmap=cm.jet,vmin=0,vmax=50)
plt.yticks([])
plt.xticks([])
plt.title('Women of Child Bearing Age 2015')
plt.colorbar()
plt.savefig(wdfigs+'pop_WOCBA',dpi=700)

A04perWOCBA=A0004/WOCBA
A04perWOCBA=A04perWOCBA[~imageMask.all(axis=1)]
A04perWOCBA=A04perWOCBA[:,~imageMask.all(axis=0)]

plt.clf()
plt.imshow(A04perWOCBA,cmap=cm.jet)
plt.yticks([])
plt.xticks([])
plt.title('Children 0-4 per Women of Child Bearing Age 2015')
plt.colorbar()
plt.savefig(wdfigs+'pop_A04perWOCBA',dpi=700)
