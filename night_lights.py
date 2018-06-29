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

def checkIfInPolygon(lon, lat,polygon):
	point = Point(lon, lat)
	return polygon.contains(point)

wddata='/Users/lilllianpetersen/iiasa/data/'
wdfigs='/Users/lilllianpetersen/iiasa/figs/'
wdvars='/Users/lilllianpetersen/iiasa/saved_vars/'

tifnlight = TIFF.open(wddata+'night_light_ethiopia/2017/2013/F182013.avg_lights_x_pct.tif')
ds=gdal.Open(wddata+'night_light_ethiopia/2017/2013/F182013.avg_lights_x_pct.tif')
width = ds.RasterXSize
height = ds.RasterYSize
gt = ds.GetGeoTransform()
minx = gt[0]
miny = gt[3] + width*gt[4] + height*gt[5] 
maxx = gt[0] + width*gt[1] + height*gt[2]
maxy = gt[3] 

pixelsize=0.004491576420598
ethMask=np.ones(shape=(width,height))
lat=np.ones(shape=(height))
lon=np.ones(shape=(width))
for w in range(width):
	lon[w]=minx+w*pixelsize
for h in range(height):
	lat[h]=miny+h*pixelsize

# read your shapefile
plt.clf()
r = shapefile.Reader(wddata+'night_light_ethiopia/2017/2013/Ethiopia_regions_2014/Ethiopia_regions_2014.shp')
rCountry = shapefile.Reader(wddata+'night_light_ethiopia/2017/2013/Ethiopia_regions_2014/ETH_outline.shp')

# Read in Regions
for i in range(len(r.shapes())):
	vars()['shapes'+str(i)] = r.shapes()[i]
	vars()['shape_points'+str(i)]=vars()['shapes'+str(i)].points
	vars()['polygon'+str(i)]=patches.Polygon(vars()['shape_points'+str(i)])
	vars()['points'+str(i)]=np.array(vars()['shape_points'+str(i)])

	plt.plot(vars()['points'+str(i)][:,0],vars()['points'+str(i)][:,1],'.')
	plt.savefig(wdfigs+'regions',dpi=600)

# Entire Country
shapeC=rCountry.shapes()[0]
shape_pointsC=shapeC.points
polygonC=patches.Polygon(shape_pointsC)
pointsC=np.array(shape_pointsC)

plt.clf()
plt.plot(pointsC[:,0],pointsC[:,1])
plt.savefig(wdfigs+'country',dpi=600)

############## Make mask for which regions the night lights are in ###############
prevRegion=0
regionMask=np.ones(shape=(len(r.shapes()),width,height),dtype=bool)
for w in range(width):
	print w
	for h in range(height):
		if checkIfInPolygon(lon[w],lat[h],polygonC)[0]==False:
			continue
		if checkIfInPolygon(lon[w],lat[h],vars()['polygon'+str(prevRegion)])[0]==True:
			regionMask[prevRegion,w,h]=False  # False=good
		else:
			for i in range(len(r.shapes())):
				if checkIfInPolygon(lon[w],lat[h],vars()['polygon'+str(i)])[0]==True:
					regionMask[i,w,h]=False  # False=good
					prevRegion=i
					break

np.save(wdvars+'regionMaskforNightLights',regionMask)
##################################################################################

#plt.clf()
#map = Basemap(llcrnrlon=minx-10,llcrnrlat=miny-10,urcrnrlon=maxx+10,urcrnrlat=maxy+10, projection='lcc',lon_0=40.49,lat_0=9.145,resolution='c')
#map.drawcoastlines(linewidth=2)
#map.drawcountries()
#map.drawmapboundary(fill_color='aqua')
#map.fillcontinents(color='#cc9966',lake_color='#99ffff')
#
#map.readshapefile(wddata+'night_light_ethiopia/2017/2013/Ethiopia_regions_2014/Ethiopia_regions_2014', name='regions', drawbounds=True, color='k',linewidth=.7)
#ax = plt.gca()
#plt.title('Ethiopia Regions')
#plt.yticks([])
#plt.xticks([])
#plt.savefig(wdfigs+'ethiopian_regions',dpi=700)

nlight = tifnlight.read_image()


nlight=np.ma.masked_array(nlight,imageMask)
nlight=nlight[~imageMask.all(axis=1)]
nlight=nlight[:,~imageMask.all(axis=0)]

plt.clf()
plt.imshow(nlight,cmap=cm.gist_gray)
plt.yticks([])
plt.xticks([])
plt.title('Night Lights Ethiopia 2013')
plt.colorbar()
plt.savefig(wdfigs+'night_lights_ethiopia_2013',dpi=700)
	
