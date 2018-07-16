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
from math import sin, cos, sqrt, atan2, radians
from geopy.geocoders import Nominatim
geolocator = Nominatim()
import geopy.distance
from scipy import ndimage
from scipy.signal import convolve2d
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

def findNearestCities(cityLonLat,gridMid,imageMask1):
	closestDist=10000*np.ones(shape=(len(lonM),len(latM)))
	closestDistDeg=100*np.ones(shape=(len(lonM),len(latM),3))
	iclosestDist=np.zeros(shape=(len(lonM),len(latM),3))
	imageMask4=np.swapaxes(imageMask3,0,1)
	R = 6373.0 #earth's radius
	for ilon in range(len(gridMid[:,0,0])):
		for ilat in range(len(gridMid[0,:,1])):
			if imageMask1[ilon,ilat]==True:
				continue
			lon,lat=gridMid[ilon,0,0],gridMid[0,ilat,1]
			for city in range(len(cityLonLat[:,0])):
				clon,clat=cityLonLat[city,0],cityLonLat[city,1]
				dist=sqrt((clon-lon)**2+(clat-lat)**2)
				if dist<closestDistDeg[ilon,ilat,0]:
					closestDistDeg[ilon,ilat,2]=closestDistDeg[ilon,ilat,1]
					iclosestDist[ilon,ilat,2]=iclosestDist[ilon,ilat,1]
					closestDistDeg[ilon,ilat,1]=closestDistDeg[ilon,ilat,0]
					iclosestDist[ilon,ilat,1]=iclosestDist[ilon,ilat,0]
					closestDistDeg[ilon,ilat,0]=dist
					iclosestDist[ilon,ilat,0]=city
					clonBest,clatBest=clon,clat
				elif dist<closestDistDeg[ilon,ilat,1]:
					closestDistDeg[ilon,ilat,2]=closestDistDeg[ilon,ilat,1]
					iclosestDist[ilon,ilat,2]=iclosestDist[ilon,ilat,1]
					closestDistDeg[ilon,ilat,1]=dist
					iclosestDist[ilon,ilat,1]=dist
				elif dist<closestDistDeg[ilon,ilat,2]:
					closestDistDeg[ilon,ilat,2]=dist
					iclosestDist[ilon,ilat,2]=dist
			
			coord1=[lat,lon]
			coord2=[clatBest,clonBest]
			closestDist[ilon,ilat]=geopy.distance.distance(coord1,coord2).km
		print np.round(100*ilon/float(len(gridMid[:,0,0])),2),'%'
	closestDistM=np.ma.masked_array(closestDist,imageMask1)
	closestDistDegM=np.ma.masked_array(closestDistDeg,imageMask4)
	iclosestDistM=np.ma.masked_array(iclosestDist,imageMask4)
	closestDistM=np.swapaxes(closestDistM,0,1)
	closestDistDegM=np.swapaxes(closestDistDegM,0,1)
	iclosestDistM=np.swapaxes(iclosestDistM,0,1)
	return closestDistM,iclosestDistM


DTYPEf = np.float64
#ctypedef np.float64_t DTYPEf_t

DTYPEi = np.int32
#ctypedef np.int32_t DTYPEi_t

#@cython.boundscheck(False) # turn of bounds-checking for entire function
#@cython.wraparound(False) # turn of bounds-checking for entire function
def replace_nans(array, max_iter, tol,kernel_size=1,method='localmean'):
    """Replace NaN elements in an array using an iterative image inpainting algorithm.
The algorithm is the following:
1) For each element in the input array, replace it by a weighted average
of the neighbouring elements which are not NaN themselves. The weights depends
of the method type. If ``method=localmean`` weight are equal to 1/( (2*kernel_size+1)**2 -1 )
2) Several iterations are needed if there are adjacent NaN elements.
If this is the case, information is "spread" from the edges of the missing
regions iteratively, until the variation is below a certain threshold.
Parameters
----------
array : 2d np.ndarray
an array containing NaN elements that have to be replaced
max_iter : int
the number of iterations
kernel_size : int
the size of the kernel, default is 1
method : str
the method used to replace invalid values. Valid options are
`localmean`, 'idw'.
Returns
-------
filled : 2d np.ndarray
a copy of the input array, where NaN elements have been replaced.
"""
    
    filled = np.empty( [array.shape[0], array.shape[1]], dtype=DTYPEf)
    kernel = np.empty( (2*kernel_size+1, 2*kernel_size+1), dtype=DTYPEf )
    
    # indices where array is NaN
    inans, jnans = np.nonzero( np.isnan(array) )
    
    # number of NaN elements
    n_nans = len(inans)
    
    # arrays which contain replaced values to check for convergence
    replaced_new = np.zeros( n_nans, dtype=DTYPEf)
    replaced_old = np.zeros( n_nans, dtype=DTYPEf)
    
    # depending on kernel type, fill kernel array
    if method == 'localmean':
      
        print 'kernel_size', kernel_size       
        for i in range(2*kernel_size+1):
            for j in range(2*kernel_size+1):
                kernel[i,j] = 1
        print kernel, 'kernel'

    elif method == 'idw':
        kernel = np.array([[0, 0.5, 0.5, 0.5,0],
                  [0.5,0.75,0.75,0.75,0.5], 
                  [0.5,0.75,1,0.75,0.5],
                  [0.5,0.75,0.75,0.5,1],
                  [0, 0.5, 0.5 ,0.5 ,0]])
        print kernel, 'kernel'      
    else:
        raise ValueError( 'method not valid. Should be one of `localmean`.')
    
    # fill new array with input elements
    for i in range(array.shape[0]):
        for j in range(array.shape[1]):
            filled[i,j] = array[i,j]

    # make several passes
    # until we reach convergence
    for it in range(max_iter):
        print 'iteration', it
        # for each NaN element
        for k in range(n_nans):
            i = inans[k]
            j = jnans[k]
            
            # initialize to zero
            filled[i,j] = 0.0
            n = 0
            
            # loop over the kernel
            for I in range(2*kernel_size+1):
                for J in range(2*kernel_size+1):
                   
                    # if we are not out of the boundaries
                    if i+I-kernel_size < array.shape[0] and i+I-kernel_size >= 0:
                        if j+J-kernel_size < array.shape[1] and j+J-kernel_size >= 0:
                                                
                            # if the neighbour element is not NaN itself.
                            if filled[i+I-kernel_size, j+J-kernel_size] == filled[i+I-kernel_size, j+J-kernel_size] :
                                
                                # do not sum itself
                                if I-kernel_size != 0 and J-kernel_size != 0:
                                    
                                    # convolve kernel with original array
                                    filled[i,j] = filled[i,j] + filled[i+I-kernel_size, j+J-kernel_size]*kernel[I, J]
                                    n = n + 1*kernel[I,J]
                                    #print n

            # divide value by effective number of added elements
            if n != 0:
                filled[i,j] = filled[i,j] / n
                replaced_new[k] = filled[i,j]
            else:
                filled[i,j] = np.nan
                
        # check if mean square difference between values of replaced
        #elements is below a certain tolerance
        print 'tolerance', np.mean( (replaced_new-replaced_old)**2 )
        if np.mean( (replaced_new-replaced_old)**2 ) < tol:
            break
        else:
            for l in range(n_nans):
                replaced_old[l] = replaced_new[l]
    
    return filled


def sincinterp(image, x,  y, kernel_size=3 ):
    """Re-sample an image at intermediate positions between pixels.
	This function uses a cardinal interpolation formula which limits the loss of information in the resampling process. It uses a limited number of neighbouring pixels.
	The new image :math:`im^+` at fractional locations :math:`x` and :math:`y` is computed as:
	.. math::
	im^+(x,y) = \sum_{i=-\mathtt{kernel\_size}}^{i=\mathtt{kernel\_size}} \sum_{j=-\mathtt{kernel\_size}}^{j=\mathtt{kernel\_size}} \mathtt{image}(i,j) sin[\pi(i-\mathtt{x})] sin[\pi(j-\mathtt{y})] / \pi(i-\mathtt{x}) / \pi(j-\mathtt{y})
	Parameters
	----------
	image : np.ndarray, dtype np.int32 the image array.
	x : two dimensions np.ndarray of floats an array containing fractional pixel row positions at which to interpolate the image
	y : two dimensions np.ndarray of floats an array containing fractional pixel column positions at which to interpolate the image
	kernel_size : int
	interpolation is performed over a ``(2*kernel_size+1)*(2*kernel_size+1)`` submatrix in the neighbourhood of each interpolation point.
	Returns
	-------
	im : np.ndarray, dtype np.float64
	the interpolated value of ``image`` at the points specified
	by ``x`` and ``y``
	"""
    
    # the output array
    r = np.zeros( [x.shape[0], x.shape[1]], dtype=DTYPEf)
          
    # fast pi
    pi = 3.1419
        
    # for each point of the output array
    for I in range(x.shape[0]):
        for J in range(x.shape[1]):
            
            #loop over all neighbouring grid points
            for i in range( int(x[I,J])-kernel_size, int(x[I,J])+kernel_size+1 ):
                for j in range( int(y[I,J])-kernel_size, int(y[I,J])+kernel_size+1 ):
                    # check that we are in the boundaries
                    if i >= 0 and i <= image.shape[0] and j >= 0 and j <= image.shape[1]:
                        if (i-x[I,J]) == 0.0 and (j-y[I,J]) == 0.0:
                            r[I,J] = r[I,J] + image[i,j]
                        elif (i-x[I,J]) == 0.0:
                            r[I,J] = r[I,J] + image[i,j] * np.sin( pi*(j-y[I,J]) )/( pi*(j-y[I,J]) )
                        elif (j-y[I,J]) == 0.0:
                            r[I,J] = r[I,J] + image[i,j] * np.sin( pi*(i-x[I,J]) )/( pi*(i-x[I,J]) )
                        else:
                            r[I,J] = r[I,J] + image[i,j] * np.sin( pi*(i-x[I,J]) )*np.sin( pi*(j-y[I,J]) )/( pi*pi*(i-x[I,J])*(j-y[I,J]))
    return r	

def moving_average_2d(data, window):
    """Moving average on two-dimensional data.
    """
    # Makes sure that the window function is normalized.
    window /= window.sum()
    # Makes sure data array is a numpy array or masked array.
    if type(data).__name__ not in ['ndarray', 'MaskedArray']:
        data = numpy.asarray(data)

    # The output array has the same dimensions as the input data 
    # (mode='same') and symmetrical boundary conditions are assumed
    # (boundary='symm').
    return convolve2d(data, window, mode='same', boundary='symm')

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
plt.imshow(pov125,cmap=cm.jet_r,vmin=0,vmax=100)
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
roads0=roadDensity>0
roads0=1*roads0
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
	distToRoads=np.load(wdvars+'nigeria_distToRoadsM.npy')
	distToRoadsM=np.ma.masked_array(distToRoads,imageMask3)
except:
	distToRoads=1-roads0 # road = 0
	distToRoads[distToRoads==1]=9999.
	
	for iroad in range(1,3):
		for i in range(300):
			print np.round(100*i/float(300),2),'%'
			for ilon in range(len(lonM)):
				for ilat in range(len(latM)):
					if imageMask1[ilon,ilat]==False:
						try:
							if distToRoads[ilon,ilat,iroad]==9999. and np.amin(distToRoads[ilon-1:ilon+2, ilat-1:ilat+2, iroad])==i:
								distToRoads[ilon,ilat,iroad]=i+1
						except:
							continue
	
	distToRoadstmp=np.swapaxes(distToRoads[:,:,:],0,1)
	distToRoadsM=np.ma.masked_array(distToRoadstmp,imageMask3)
	whereBad=np.where(distToRoadsM[:,:]==9999)
	for i in range(whereBad[0].shape[0]):
		ilat=whereBad[0][i]
		ilon=whereBad[1][i]
		iroad=whereBad[2][i]
		if ilat==0:
			distToRoadsM[ilat,ilon,iroad]=np.amin(distToRoadsM[ilat:ilat+2,ilon-1:ilon+2,iroad])
		if ilon==0:
			distToRoadsM[ilat,ilon,iroad]=np.amin(distToRoadsM[ilat-1:ilat+2,ilon:ilon+2,iroad])
	np.save(wdvars+'nigeria_distToRoadsM',np.array(distToRoadsM))

plt.clf()
plt.imshow(np.clip(distToRoadsM[:,:,1],0,15),cmap=cm.jet_r)
plt.title('e^ Dist To Primary Roads')
plt.colorbar()
plt.savefig(wdfigs+'dist_to_primary_roads_clipped',dpi=700)

plt.clf()
plt.imshow(np.exp2(distToRoadsM[:,:,1]),cmap=cm.jet_r,vmin=0,vmax=10)
plt.title('e^ Dist To Primary Roads')
plt.colorbar()
plt.savefig(wdfigs+'dist_to_primary_roads',dpi=700)

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

#pov125M=np.ma.masked_array(pov125,np.swapaxes(distToRoadsMask[:,:,0],0,1))
#var=['Primary','Secondary','Tertiary']
#Corr=np.zeros(shape=(6))
#slope=np.zeros(shape=(6))
#bInt=np.zeros(shape=(6))
#xMulti=np.zeros(shape=(1079449,3))
#for i in range(3):
#	x=np.ma.compressed(distToRoadsM[:,:,i])
#	xMulti[:,i]=x
#	ydata=np.ma.compressed(pov125M)
#	yMulti=ydata
#	Corr[i]=corr(x,ydata)
#	slope[i],bInt[i]=np.polyfit(x,ydata,1)
#	yfit=slope[i]*x+bInt[i]
#	plt.clf()
#	plt.figure(21,figsize=(8,6))
#	plt.plot(x,ydata,'*',markersize=1)
#	plt.plot(x,yfit,'g-')
#	plt.xlabel('Dist to '+var[i]+' Roads')
#	plt.ylabel('% Below 1.25$/day')
#	plt.title('Distance to '+var[i]+' Roads and Poverty %, Corr = '+str(np.round(Corr[i],2)))
#	plt.savefig(wdfigs+'roaddist_'+var[i]+'_pov_corr',dpi=700)
#
#roadDensityM2=np.array(roadDensityM2)
#roadDensityM2mask=roadDensityM==0
#roadDensityM2=np.ma.masked_array(roadDensityM,roadDensityM2mask)
#for i in range(3):
#	pov125M=np.ma.masked_array(pov125,roadDensityM2mask[:,:,i])
#	x=np.ma.compressed(roadDensityM2[:,:,i])
#	print x.shape
#	ydata=np.ma.compressed(pov125M)
#	Corr[i+3]=corr(x,ydata)
#	slope[i+3],bInt[i+3]=np.polyfit(x,ydata,1)
#	yfit=slope[i+3]*x+bInt[i+3]
#	plt.clf()
#	plt.figure(21,figsize=(8,6))
#	plt.plot(x,ydata,'*',markersize=1)
#	plt.plot(x,yfit,'g-')
#	plt.xlabel('Density of '+var[i]+' Roads')
#	plt.ylabel('% Below 1.25$/day')
#	plt.title('Density of '+var[i]+' Roads and Poverty %, Corr = '+str(np.round(Corr[i+3],2)))
#	plt.savefig(wdfigs+'roaddensity_'+var[i]+'_pov_corr',dpi=700)

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
#from sklearn import linear_model
#
#clf=linear_model.LinearRegression()
#clf.fit(xMulti,yMulti)
#exit()

#np.save(wdvars+'Illinois/xMulti',xMulti)
#np.save(wdvars+'Illinois/ydataMulti',ydata)
#
#
################ Claculate Distance to Cities ###############

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
	#closestDistDegM=np.load(wdvars+'closestDistDegM_to_nigerian_cities.npy')
except:
	gridMid=createMidpointGrid(grid,pixelsize)
	closestDistM,iclosestDistM=findNearestCities(cityLonLat,gridMid,imageMask1)
	np.save(wdvars+'closestDist_to_nigerian_cities',np.array(closestDistM))
	np.save(wdvars+'iclosestDist_to_nigerian_cities',np.array(iclosestDistM))

plt.clf()
plt.imshow(closestDistM,cmap=cm.viridis_r)
plt.colorbar()
plt.title('Distance to Nearest City')
plt.savefig(wdfigs+'closestDist',dpi=700)

popClosestDist=np.zeros(shape=(iclosestDistM.shape))
popClosestDistScaled=np.zeros(shape=(iclosestDistM.shape))
for ilat in range(len(iclosestDistM[:,0])):
	for ilon in range(len(iclosestDistM[0,:])):
		if imageMask2[ilat,ilon]==False:
			#popClosestDistScaled[ilat,ilon]=popScaled[int(iclosestDistM[ilat,ilon])]
			popClosestDist[ilat,ilon]=cityPop[int(iclosestDistM[ilat,ilon])]
popClosestDistM=np.ma.masked_array(popClosestDist,imageMask2)

######## Index try 1 ##########
#dmin=0
#dmax=np.amax(cityPop)
#dif=dmax-dmin
#popScaled=((cityPop-dmin)/dif)
#popScaled=1-popScaled
#
#
#dmin=0
#dmax=np.amax(closestDistM)
#dif=dmax-dmin
#max1=1
#min1=.3
#closestDistScaled=((closestDistM-dmin)/dif)*(max1-min1)+min1
#distToCityIndex=closestDistScaled*popClosestDistScaled


'''
############## Playing around with data ##############
# Smooth the data
indexL=np.log(index)
indexLS=indexL
indexLS[imageMask2==True]=6
win = np.ones((30, 30))
indexLS=moving_average_2d(indexL,win)
indexLS=np.ma.masked_array(indexLS,imageMask2)
x=np.ma.compressed(indexLS)
ydata=np.ma.compressed(pov125)
Corr=corr(x,ydata)
slope,bInt=np.polyfit(x,ydata,1)
yfit=slope*x+bInt
plt.clf()
plt.plot(x,ydata,'*',markersize=.2)
plt.plot(x,yfit)
plt.title('Index and Pov125, Corr = '+str(round(Corr,2)))
plt.ylabel('pov125')
plt.xlabel('population/distance to city')
plt.savefig(wdfigs+'popOverDistanceAndPovCorr',dpi=700)

iclosestDistM[imageMask2==True]=-99
iclosestDistM=np.array(iclosestDistM)
mask=iclosestDistM==0
mask=1-mask
mask=np.array(mask,dtype=bool)
roadDensityMtmp=np.ma.masked_array(roadDensityM[:,:,2],mask)
roadDensityMtmp=roadDensityMtmp[~mask.all(axis=1)]
roadDensityMtmp=roadDensityMtmp[:,~mask.all(axis=0)]
plt.clf()
plt.imshow(roadDensityMtmp)
plt.colorbar()
plt.savefig(wdfigs+'example_roads_teritiary',dpi=700)

pov125tmp=np.ma.masked_array(pov125,mask)
pov125tmp=pov125tmp[~mask.all(axis=1)]
pov125tmp=pov125tmp[:,~mask.all(axis=0)]
plt.clf()
plt.imshow(pov125tmp)
plt.colorbar()
plt.savefig(wdfigs+'example_pov',dpi=700)

closestDisttmp=np.ma.masked_array(closestDistM,mask)
closestDisttmp=closestDisttmp[~mask.all(axis=1)]
closestDisttmp=closestDisttmp[:,~mask.all(axis=0)]
plt.clf()
plt.imshow(closestDisttmp)
plt.colorbar()
plt.savefig(wdfigs+'example_dist',dpi=700)

distToRoadsMtmp=np.ma.masked_array(distToRoadsM[:,:,0],mask)
distToRoadsMtmp=distToRoadsMtmp[~mask.all(axis=1)]
distToRoadsMtmp=distToRoadsMtmp[:,~mask.all(axis=0)]
plt.clf()
plt.imshow(distToRoadsMtmp)
plt.colorbar()
plt.savefig(wdfigs+'example_dist_roads',dpi=700)

distToRoadsS=(distToRoadsMtmp-np.amin(distToRoadsMtmp))/(np.amax(distToRoadsMtmp)-np.amin(distToRoadsMtmp))
roadDensityS=(roadDensityMtmp-np.amin(roadDensityMtmp))/(np.amax(roadDensityMtmp)-np.amin(roadDensityMtmp))*.7
closestDistS=(closestDisttmp-np.amin(closestDisttmp))/(np.amax(closestDisttmp)-np.amin(closestDisttmp))*.7
index1=distToRoadsS-roadDensityS+closestDistS
plt.clf()
plt.imshow(index1)
plt.colorbar()
plt.title('DistToRoadsP - TRoadDensity + DistToCity')
plt.savefig(wdfigs+'example_dist-density',dpi=700)
print corr(np.ma.compressed(index1),np.ma.compressed(pov125tmp))
x=np.ma.compressed(index1)
ydata=np.ma.compressed(pov125tmp)
Corr=corr(x,ydata)
slope,bInt=np.polyfit(x,ydata,1)
yfit=slope*x+bInt
plt.clf()
plt.plot(x,ydata,'*')
plt.plot(x,yfit)
plt.title('Index1 and Pov125, Corr = '+str(round(Corr,2)))
plt.ylabel('pov125')
plt.xlabel('DistToRoadsP - TRoadDensity + DistToCity')
plt.savefig(wdfigs+'index1_and_pov125',dpi=700)
'''

###########################################
# Night Lights
###########################################
stableTif=TIFF.open(wddata+'nigeria_lights/nigeriaStable.tif',mode='r')
avgVisTif=TIFF.open(wddata+'nigeria_lights/nigeriaAvgVis.tif',mode='r')

ds=gdal.Open(wddata+'nigeria_lights/nigeriaAvgVis.tif')
width = ds.RasterXSize
height = ds.RasterYSize
gt = ds.GetGeoTransform()
minx = gt[0]
miny = gt[3] + width*gt[4] + height*gt[5] 
maxx = gt[0] + width*gt[1] + height*gt[2]
maxy = gt[3] 
pixelsizel=abs(gt[-1])

latl=np.ones(shape=(height))
lonl=np.ones(shape=(width))
for w in range(width):
	lonl[w]=minx+w*pixelsizel
for h in range(height):
	latl[h]=miny+h*pixelsizel

stableSmall=stableTif.read_image()
avgVisSmall=avgVisTif.read_image()

latl=np.radians(latl)
lonl=np.radians(lonl)
lutAvgVis=RectSphereBivariateSpline(latl, lonl, avgVisSmall)
lutStable=RectSphereBivariateSpline(latl, lonl, stableSmall)

newLats=np.radians(latM[::-1])
newLons=np.radians(lonM)
newLats,newLons=np.meshgrid(newLats,newLons)
avgVisO=lutAvgVis.ev(newLats.ravel(),newLons.ravel()).reshape((len(lonM),len(latM))).T
stable=lutStable.ev(newLats.ravel(),newLons.ravel()).reshape((len(lonM),len(latM))).T

avgVis=np.log(avgVisO)
avgVis[np.isnan(avgVis)]=0

avgVis=(avgVis-np.amin(avgVis))/(np.amax(avgVis)-np.amin(avgVis))
stable=(stable-np.amin(stable))/(np.amax(stable)-np.amin(stable))
avgVis=np.clip(avgVis,0.94,1)
avgVis=(avgVis-np.amin(avgVis))/(np.amax(avgVis)-np.amin(avgVis))

avgVisS=moving_average_2d(avgVis,np.ones((5, 5)))
avgVisS=np.ma.masked_array(avgVisS,imageMask2)
avgVisS=(avgVisS-np.amin(avgVisS))/(np.amax(avgVisS)-np.amin(avgVisS))

plt.clf()
plt.imshow(avgVisS,cmap=cm.jet)
plt.colorbar()
plt.title('Avg Visual Lights Smoothed')
plt.savefig(wdfigs+'nigeria_avgVisS',dpi=700)

plt.clf()
plt.imshow(stable,cmap=cm.jet)
plt.colorbar()
plt.title('Stable Lights')
plt.savefig(wdfigs+'nigeria_stable_lights',dpi=700)

#plt.clf()
#plt.imshow(avgVis,cmap=cm.jet)
#plt.colorbar()
#plt.title('Avg Visual Lights')
#plt.savefig(wdfigs+'nigeria_avgVis',dpi=700)
#
#avgVis=np.ma.masked_array(avgVis,imageMask2)
stable=np.ma.masked_array(stable,imageMask2)


###########################################
# Make Index
###########################################
#### For all of Nigeria ####
#cityDensity=closestDistDegM[:,:,0]+closestDistDegM[:,:,1]+closestDistDegM[:,:,2]
#cityDensity=(cityDensity-np.amin(cityDensity))/(np.amax(cityDensity)-np.amin(cityDensity))
#plt.clf()
#plt.imshow(cityDensity,cmap=cm.hot_r)
#plt.title('city density')
#plt.colorbar()
#plt.savefig(wdfigs+'city_density',dpi=700)
#
#cityDensityPop=np.zeros(shape=(len(latM),len(lonM)))
#for ilat in range(len(latM)):
#	for ilon in range(len(lonM)):
#		if imageMask2[ilat,ilon]==False:
#			cityDensityPop[ilat,ilon]=cityPop[int(iclosestDistM[ilat,ilon,0])]+cityPop[int(iclosestDistM[ilat,ilon,1])]+cityPop[int(iclosestDistM[ilat,ilon,2])]
#cityDensityPop=cityDensityPop/cityDensity
#cityDensityPop=np.log(cityDensityPop)
#cityDensityPopM=np.ma.masked_array(cityDensityPop,imageMask2)
#plt.clf()
#plt.imshow(cityDensityPopM,cmap=cm.hot_r)
#plt.title('City Density by Population')
#plt.colorbar()
#plt.savefig(wdfigs+'city_density_by_pop',dpi=700)

#distToPRoadsS=(distToRoadsM[:,:,0]-np.amin(distToRoadsM[:,:,0]))/(np.amax(distToRoadsM[:,:,0])-np.amin(distToRoadsM[:,:,0]))
#distToSRoadsS=(distToRoadsM[:,:,1]-np.amin(distToRoadsM[:,:,1]))/(np.amax(distToRoadsM[:,:,1])-np.amin(distToRoadsM[:,:,1]))
# Primary
distToRoadsClipped=np.clip(distToRoadsM[:,:,0],0,15)
distToPRoadsS=(distToRoadsClipped[:,:]-np.amin(distToRoadsClipped[:,:]))/(np.amax(distToRoadsClipped[:,:])-np.amin(distToRoadsClipped[:,:]))
# Secondary
distToRoadsClipped=np.clip(distToRoadsM[:,:,1],0,10)
distToSRoadsS=(distToRoadsClipped[:,:]-np.amin(distToRoadsClipped[:,:]))/(np.amax(distToRoadsClipped[:,:])-np.amin(distToRoadsClipped[:,:]))
# Teritiary
distToRoadsClipped=np.clip(distToRoadsM[:,:,2],0,3)
distToTRoadsS=(distToRoadsClipped[:,:]-np.amin(distToRoadsClipped[:,:]))/(np.amax(distToRoadsClipped[:,:])-np.amin(distToRoadsClipped[:,:]))

TroadDensityS=(roadDensityM[:,:,2]-np.amin(roadDensityM[:,:,2]))/(np.amax(roadDensityM[:,:,2])-np.amin(roadDensityM[:,:,2]))
SroadDensityS=(roadDensityM[:,:,1]-np.amin(roadDensityM[:,:,1]))/(np.amax(roadDensityM[:,:,1])-np.amin(roadDensityM[:,:,1]))
#roadDensityS=popClosestDistM/roadDensityS
closestDistS=(closestDistM-0)/(np.amax(closestDistM)-0)
#closestDistS=popClosestDistM/closestDistS
try:
	popClosestDistS=np.load(wdvars+'nigeria_popClosestDistSmoothed.npy')
	popClosestDistS=np.ma.masked_array(popClosestDistS,imageMask2)
except:
	popClosestDist1=popClosestDistM.filled(np.NaN)
	popClosestDistFilled = replace_nans(popClosestDist1, 10, 0.5, 2, method='localmean')
	popClosestDistS=np.log(popClosestDistFilled)
	popClosestDistS=moving_average_2d(popClosestDistS,np.ones((40, 40)))
	popClosestDistS=np.ma.masked_array(popClosestDistS,imageMask2)
	popClosestDistSmoothed=popClosestDistS
	popClosestDistSmoothed[imageMask2==True]=0
	popClosestDistSmoothed=np.array(popClosestDistSmoothed)
	np.save(wdvars+'nigeria_popClosestDistSmoothed',popClosestDistSmoothed)
dmin=np.amin(popClosestDistS)
popClosestDistSS=(popClosestDistS-dmin)/(np.amax(popClosestDistS)-dmin)+0.01

index=popClosestDistM/closestDistM
index=np.log(index)
index=(index-np.amin(index))/(np.amax(index)-np.amin(index))
index=1-index
plt.clf()
plt.imshow(index,cmap=cm.jet_r)
plt.title('Population of Closest City / Distance')
plt.colorbar()
plt.savefig(wdfigs+'index',dpi=700)


TroadDensityS1=TroadDensityS.filled(np.NaN)
TroadDensityFilled=replace_nans(TroadDensityS1,2,0.5, 2, method='localmean')
TroadDensity=moving_average_2d(TroadDensityFilled,np.ones((3,3)))
TroadDensity=np.ma.masked_array(TroadDensity,imageMask2)
TroadDensity=(TroadDensity-np.amin(TroadDensity))/(np.amax(TroadDensity)-np.amin(TroadDensity))

i1=1.5*distToPRoadsS+distToSRoadsS+distToTRoadsS
i2=-TroadDensityS-SroadDensityS
i1=(i1-np.amin(i1))/(np.amax(i1)-np.amin(i1))
i2=(i2-np.amin(i2))/(np.amax(i2)-np.amin(i2))
indext=i1+i2
indext=(indext-np.amin(indext))/(np.amax(indext)-np.amin(indext))
#indext=-roadDensityS
plt.clf()
plt.imshow(indext,cmap=cm.hot_r)
plt.colorbar()
plt.title('DistToRoadsP - TRoadDensity')
plt.savefig(wdfigs+'indext',dpi=700)

indext1=np.log(popClosestDistSS/closestDistS)
indext1=(indext1-np.amin(indext1))/(np.amax(indext1)-np.amin(indext1))
indext1=1-indext1
#indext1=closestDistS
plt.clf()
plt.imshow(indext1,cmap=cm.hot_r)
plt.colorbar()
plt.title('popCity/DistToCity')
plt.savefig(wdfigs+'indext1',dpi=700)

#index1=distToPRoadsS-roadDensityS+closestDistS-popClosestDistSS
#index1=2*indext+indext1
index1=1.25*index+indext-3*avgVisS
index1=(index1-np.amin(index1))/(np.amax(index1)-np.amin(index1))
#index11=index1.filled(np.NaN)
#index1Filled=replace_nans(index11,3,0.5, 2, method='localmean')
#index1s=moving_average_2d(index1Filled,np.ones((5,5)))
#index1s=np.ma.masked_array(index1s,imageMask2)
#index1s=(index1s-np.amin(index1s))/(np.amax(index1s)-np.amin(index1s))

###########################################
# Scale and find Corr
###########################################
scaleIndex=10.

height,width=int(np.round(len(latM)/scaleIndex)),int(np.round(len(lonM)/scaleIndex))
newLats,newLons=np.ones(shape=(height)),np.ones(shape=(width))
pixelsizeNew=pixelsize*scaleIndex
latMr=latM[::-1]
for w in range(width):
	newLons[w]=lonM[0]+w*pixelsizeNew
for h in range(height):
	newLats[h]=latMr[0]+h*pixelsizeNew

gridNew = np.zeros(shape=(len(newLons),len(newLats),2))
for x in range(len(newLons)):
	gridNew[x,:,0]=newLons[x]
for y in range(len(newLats)):
	gridNew[:,y,1]=newLats[y]
grid1=grid[:,::-1,:]

## Plot old and new grid
plt.clf()
m=int(scaleIndex*3)
plt.plot(np.ma.compressed(grid1[:m,:m,0]),np.ma.compressed(grid1[:m,:m,1]),'*k')
plt.plot(np.ma.compressed(gridNew[:3,:3,0]),np.ma.compressed(gridNew[:3,:3,1]),'*r')
plt.savefig(wdfigs+'old_new_grid.pdf')
##

oldLat,oldLon=np.radians(latM[::-1]),np.radians(lonM)
lutIndex1=RectSphereBivariateSpline(oldLat, oldLon, index1)
lutPov=RectSphereBivariateSpline(oldLat, oldLon, pov125)
lutMask=RectSphereBivariateSpline(oldLat, oldLon, imageMask2)

newLats,newLons=np.radians(newLats),np.radians(newLons)
newLats,newLons=np.meshgrid(newLats,newLons)

urbanIndex=lutIndex1.ev(newLats.ravel(),newLons.ravel()).reshape((width,height)).T
pov125Course=lutPov.ev(newLats.ravel(),newLons.ravel()).reshape((width,height)).T
maskCourse=lutMask.ev(newLats.ravel(),newLons.ravel()).reshape((width,height)).T
maskCourse=np.array(np.round(maskCourse,0),dtype=bool)

urbanIndex=np.ma.masked_array(urbanIndex,maskCourse)
urbanIndex=(urbanIndex-np.amin(urbanIndex))/(np.amax(urbanIndex)-np.amin(urbanIndex))
pov125Course=np.ma.masked_array(pov125Course,maskCourse)

plt.clf()
plt.imshow(urbanIndex,cmap=cm.jet_r)
plt.colorbar()
plt.title('Urban/Rural Index Based on Roads and Population')
plt.savefig(wdfigs+'urbanIndex',dpi=700)

plt.clf()
plt.imshow(pov125Course,cmap=cm.jet_r)
plt.colorbar()
plt.title('Percent Under $1.25/Day')
plt.savefig(wdfigs+'pov125Course',dpi=700)
print corr(np.ma.compressed(urbanIndex),np.ma.compressed(pov125Course))

x=np.ma.compressed(index1)
ydata=np.ma.compressed(pov125)
Corr=corr(x,ydata)
slope,bInt=np.polyfit(x,ydata,1)
yfit=slope*x+bInt
plt.clf()
plt.plot(x,ydata,'.',markersize=0.0002)
plt.plot(x,yfit)
plt.title('Index1 and Pov125, Corr = '+str(round(Corr,2)))
plt.ylabel('pov125')
plt.xlabel('distToPRoads+distToSRoads-(TroadDensity+SroadDensity)+log(popClosestDistM/closestDistM)')
plt.savefig(wdfigs+'index1_and_pov125_nigeria',dpi=700)
exit()


plt.clf()
plt.imshow(index,cmap=cm.hot_r,norm=colors.LogNorm())
plt.colorbar()
plt.title('Rural Index by Cities')
plt.savefig(wdfigs+'index',dpi=700)


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
