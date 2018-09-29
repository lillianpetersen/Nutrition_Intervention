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
from scipy.stats import norm, mode
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
from sklearn.ensemble import RandomForestRegressor
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

def findGridDataLays(shapePoints,gridMid):
	'''Convert Vector To Raster'''
	griddedHits=np.zeros(shape=(len(grid[:,0]),len(grid[0,:])))
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

def findNearestThree(cityLonLat,gridMid,imageMask1):
	lenLon,lenLat=len(gridMid[:,0,0]),len(gridMid[0,:,0])
	closestDist=10000*np.ones(shape=(lenLon,lenLat))
	closestDistDeg=100*np.ones(shape=(lenLon,lenLat,3))
	iclosestDist=np.zeros(shape=(lenLon,lenLat,3))
	imageMask4=np.zeros(shape=(lenLon,lenLat,3))
	for i in range(3):
		imageMask4[:,:,i]=imageMask1
	R = 6373.0 #earth's radius
	for ilon in range(lenLon):
		for ilat in range(lenLat):
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

################################################

wddata='/Users/lilllianpetersen/iiasa/data/'
wdfigs='/Users/lilllianpetersen/iiasa/figs/'
wdvars='/Users/lilllianpetersen/iiasa/saved_vars/'

makePlots=False

##############################################
# Read Poverty
##############################################
i=1 #country=Nigeria
country='Nigeria'

tif125 = TIFF.open(wddata+'poverty/'+country+'/nga10povcons125.tif',mode='r')
ds=gdal.Open(wddata+'poverty/'+country+'/nga10povcons125.tif')
width = ds.RasterXSize
height = ds.RasterYSize
gt = ds.GetGeoTransform()
minx = gt[0]
miny = gt[3] + width*gt[4] + height*gt[5] 
maxx = gt[0] + width*gt[1] + height*gt[2]
maxy = gt[3] 
pixelsize=abs(gt[-1])

# lat and lon
lat=np.ones(shape=(height))
lon=np.ones(shape=(width))
for w in range(width):
	lon[w]=minx+w*pixelsize
for h in range(height):
	lat[h]=miny+h*pixelsize
lat=lat[::-1] # reverse the order

pov125 = tif125.read_image()*100
imageMask=pov125<0

pov125=np.ma.masked_array(pov125,imageMask)
pov125=pov125[~imageMask.all(axis=1)]
pov125=pov125[:,~imageMask.all(axis=0)]

lonp=np.ma.masked_array(lon,imageMask.all(axis=0))
lonp=np.ma.compressed(lonp)
latp=np.ma.masked_array(lat,imageMask.all(axis=1))
latp=np.ma.compressed(latp)

## Masks
imageMask21km=imageMask # imageMask2=(lat,lon)
imageMask21km=imageMask21km[~imageMask21km.all(axis=1)]
imageMask21km=imageMask21km[:,~imageMask21km.all(axis=0)]
imageMask11km=np.swapaxes(imageMask21km,0,1) # imageMask1=(lon,lat)
imageMask31km=np.zeros(shape=(imageMask21km.shape[0],imageMask21km.shape[1],10),dtype=bool)
for i in range(10):
	imageMask31km[:,:,i]=imageMask21km # imageMask2=(lat,lon,3)


##############################################
# Gridded Malnutrition
##############################################
print 'Malnutrition' 
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

mal=ds.ReadAsArray()
mal[mal>1]=0.0308212 # The average of this weird pixel's neighbors

mal=mal[:1740]
latm=latm[:1740]
mal=mal[:,185:]
lonm=lonm[185:]

mal=mal[latm<latp[0]]
latm=latm[latm<latp[0]]
mal=mal[latm>latp[-1]]
latm=latm[latm>latp[-1]]

mal=mal[:,lonm>lonp[0]]
lonm=lonm[lonm>lonp[0]]
mal=mal[:,lonm<lonp[-1]]
lonm=lonm[lonm<lonp[-1]]

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

##### 5km mask #####
latl=np.radians(latp[::-1])
lonl=np.radians(lonp)
lut=RectSphereBivariateSpline(latl, lonl, imageMask21km)

newLats,newLons=np.meshgrid(np.radians(latm[::-1]),np.radians(lonm))
imageMask2=lut.ev(newLats.ravel(),newLons.ravel()).reshape((len(lonm),len(latm))).T
imageMask2=np.round(imageMask2,0)
imageMask2[mal==0]=1


imageMask3=np.zeros(shape=(imageMask2.shape[0],imageMask2.shape[1],10),dtype=bool)
for i in range(10):
	imageMask3[:,:,i]=imageMask2 # imageMask2=(lat,lon,3)

## grid
grid = np.zeros(shape=(len(latm),len(lonm),2))
for x in range(len(latm)):
	grid[x,:,0]=latm[x]
for y in range(len(lonm)):
	grid[:,y,1]=lonm[y]

gridMid=createMidpointGrid(grid,pixelsize) # lon lat


##############################################
# Gridded Population
##############################################
print 'Population' 
ds=gdal.Open(wddata+'population/gpw_v4_basic_demographic_characteristics_rev10_a000_004bt_2010_cntm_30_sec.tif')
#dsM=gdal.Open(wddata+'population/worldpop/AFR_PPP_A0004_M_2015_adj_v5.tif')
#dsF=gdal.Open(wddata+'population/worldpop/AFR_PPP_A0004_F_2015_adj_v5.tif')
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

#worldPopM=dsM.ReadAsArray()
#worldPopF=dsF.ReadAsArray()
#worldPop=worldPopM+worldPopF
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

plt.clf()
plt.imshow(pop,cmap=cm.gist_ncar_r,vmax=2000)
plt.title('gridded population')
plt.colorbar()
plt.savefig(wdfigs+'pop5km',dpi=900)

# This is the important line	
pop[pop<0]=0

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
pop5km=np.ma.masked_array(pop5km,imageMask2)

plt.clf()
plt.imshow(pop5km,cmap=cm.gist_ncar_r,vmax=2000)
plt.title('gridded population')
plt.colorbar()
plt.savefig(wdfigs+'pop5km',dpi=700)

malnumber=mal*pop5km
malnumber=np.ma.masked_array(malnumber,imageMask2)
mal=np.ma.masked_array(mal,imageMask2)

plt.clf()
plt.imshow(malnumber,cmap=cm.nipy_spectral_r,vmax=500)
plt.title('malnutrition number')
plt.colorbar()
plt.savefig(wdfigs+'malnumber',dpi=900)

plt.clf()
plt.imshow(mal,cmap=cm.jet,vmax=0.3)
plt.title('gridded wasting 2015')
plt.colorbar()
plt.savefig(wdfigs+'wasting',dpi=700)


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

travelWorld=traveltif.read_image()
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
travelbi=np.array(travel)
travelbi[travel!=0]=1
travelbi=1-travelbi

### 5km ###
latl=np.radians(latt[::-1])
lonl=np.radians(lont)
lut=RectSphereBivariateSpline(latl, lonl, travel)

newLats,newLons=np.meshgrid(np.radians(latm[::-1]),np.radians(lonm))
travel=lut.ev(newLats.ravel(),newLons.ravel()).reshape((len(lonm),len(latm))).T

travel=np.ma.masked_array(travel,imageMask2)

plt.clf()
plt.imshow(travel,cmap=cm.gist_ncar_r,vmax=200)
plt.title('Travel Time')
plt.colorbar()
plt.savefig(wdfigs+'travel',dpi=700)

######################################################
# Roads
######################################################
print 'Roads'
try:
	roadDensity=np.load(wdvars+'nigeria_raodDensity.npy')
except:
	print 'calculating roadDensity: try command failed'

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

	roadDensity=np.zeros(shape=(len(latm),len(lonm),3))
	
	for i in range(len(records1)):
		print i,np.round(100*i/float(len(records1)),1)
		shapePointsLonLat=records1[i].shape.points
		shapePointsLonLat=np.array(shapePointsLonLat)
		shapePoints=np.array(shapePointsLonLat)
		shapePoints[:,0]=shapePointsLonLat[:,1]
		shapePoints[:,1]=shapePointsLonLat[:,0]
		midpointHits,dist=findGridAndDistance(shapePoints,gridMid)
		midpointHits=np.array(midpointHits,dtype=int)
		for j in range(len(midpointHits)):
			roadDensity[midpointHits[j,0],midpointHits[j,1],0]+=dist[j]
	
	for i in range(len(records2)):
		print i,np.round(100*i/float(len(records2)),1)
		shapePointsLonLat=records2[i].shape.points
		shapePointsLonLat=np.array(shapePointsLonLat)
		shapePoints=np.array(shapePointsLonLat)
		shapePoints[:,0]=shapePointsLonLat[:,1]
		shapePoints[:,1]=shapePointsLonLat[:,0]
		midpointHits,dist=findGridAndDistance(shapePoints,gridMid)
		midpointHits=np.array(midpointHits,dtype=int)
		for j in range(len(midpointHits)):
			roadDensity[midpointHits[j,0],midpointHits[j,1],1]+=dist[j]
	
	for i in range(len(records3)):
		print i,np.round(100*i/float(len(records3)),2)
		shapePointsLonLat=records3[i].shape.points
		shapePointsLonLat=np.array(shapePointsLonLat)
		shapePoints=np.array(shapePointsLonLat)
		shapePoints[:,0]=shapePointsLonLat[:,1]
		shapePoints[:,1]=shapePointsLonLat[:,0]
		midpointHits,dist=findGridAndDistance(shapePoints,gridMid)
		midpointHits=np.array(midpointHits,dtype=int)
		for j in range(len(midpointHits)):
			roadDensity[midpointHits[j,0],midpointHits[j,1],2]+=dist[j]
	
	np.save(wdvars+'nigeria_raodDensity',roadDensity)

roadDensity=np.ma.masked_array(roadDensity,imageMask3[:,:,:3])

if makePlots:
	plt.clf()
	map = Basemap(llcrnrlon=lonm[0],llcrnrlat=np.amin(latm),urcrnrlon=lonm[-1],urcrnrlat=np.amax(latm), projection='lcc',lon_0=(lonm[-1]+lonm[0])/2,lat_0=(latm[-1]+latm[0])/2,resolution='i')
	map.drawcoastlines(linewidth=1.5)
	map.drawcountries()
	map.readshapefile(wddata+'openstreetmap/nigeria/openstreetmap/nga_trs_roads_osm', name='roads', drawbounds=True)
	ax = plt.gca()
	plt.title('Nigerian Roads from Open Street Map')
	plt.savefig(wdfigs+'openstreetmap_nigeria',dpi=700)

	roadDensity1=np.array(roadDensity[:])
	roadDensity1[roadDensity>0]=1
	roadDensity1=np.ma.masked_array(roadDensity1,roadDensity==0)
	plt.clf()
	plt.figure(1,figsize=(3,4))
	plt.imshow((roadDensity1!=0)[:,:,2],cmap=cm.YlGn)
	plt.imshow((roadDensity1!=0)[:,:,1],cmap=cm.RdPu)
	plt.imshow((roadDensity1!=0)[:,:,0],cmap=cm.binary)
	plt.colorbar()
	plt.xticks([])
	plt.yticks([])
	plt.title('Nigeria Roads')
	plt.savefig(wdfigs+'nigeria_roads',dpi=800)
	
plt.clf()
plt.imshow(roadDensity[:,:,2],cmap=cm.nipy_spectral_r)
plt.colorbar()
plt.xticks([])
plt.yticks([])
plt.title('Nigeria Teritiary Road Density')
plt.savefig(wdfigs+'nigeria_teritiary_road_density',dpi=800)

###### Dist To Roads ######
try:
	distToRoads=np.load(wdvars+'nigeria_distToRoads.npy')
	distToRoads=np.ma.masked_array(distToRoads,imageMask3[:,:,:3])
except:
	print 'calculating distToRoads: try command failed'
	roads0=roadDensity>0
	roads0=1*roads0
	
	distToRoads=1-roads0 # road = 0
	distToRoads[distToRoads==1]=9999.
	
	distToRoads=np.zeros(shape=(len(latm),len(lonm),3))
	for iroad in range(3):
		roadsMidPoints=gridMid[roads0[:,:,iroad]==1,:]
		print 'finding distance to roads',iroad
		distToRoads[:,:,iroad],iclosestRoad=findNearest(roadsMidPoints,gridMid,imageMask2)
	np.save(wdvars+'nigeria_distToRoads',np.array(distToRoads))
	distToRoads=np.ma.masked_array(distToRoads,imageMask3[:,:,:3])
		
plt.clf()
plt.imshow(np.clip(distToRoads[:,:,0],0,50),cmap=cm.jet_r)
plt.title('Dist To Primary Roads')
plt.colorbar()
plt.savefig(wdfigs+'dist_to_primary_roads',dpi=700)

#############################################
# Dist to Cities
#############################################
print 'Cities'
###### City Population ######
try:
	cityPop=np.load(wdvars+'NigeriaCityPop.npy')
	cityLatLon=np.load(wdvars+'NigeriaCityLatLon.npy')
except:
	print 'calculating cityPop: try command failed'
	f = open(wddata+'openstreetmap/nigeria/nigeria_cities_pop.csv','r')
	
	cityPop=np.zeros(shape=(400))
	cityLatLon=np.zeros(shape=(400,2))
	city=[]
	i=-1
	for line in f:
		i+=1
		tmp=line.split(',')
		citytmp=tmp[0]
		city.append(citytmp)
		location = geolocator.geocode(citytmp+', Nigeria')
		cityLatLon[i]=location.latitude,location.longitude
		cityPop[i]=float(tmp[1])
	np.save(wdvars+'NigeriaCityPop',cityPop)
	np.save(wdvars+'NigeriaCityLatLon',cityLatLon)

###### Dist To Cities and Pop of Closest City ######
try:
	closestDistCities=np.load(wdvars+'closestDist_to_nigerian_cities.npy')
	iclosestDistCities=np.load(wdvars+'iclosestDist_to_nigerian_cities.npy')
	popClosestDistCities=np.load(wdvars+'nigeria_popClosestDistCities.npy')
except:
	print 'calculating closestDistCities: try command failed'
	closestDistCities,iclosestDistCities=findNearest(cityLatLon,gridMid,imageMask2)
	np.save(wdvars+'closestDist_to_nigerian_cities',np.array(closestDistCities))
	np.save(wdvars+'iclosestDist_to_nigerian_cities',np.array(iclosestDistCities))

	plt.clf()
	plt.imshow(closestDistCities,cmap=cm.viridis_r)
	plt.colorbar()
	plt.yticks([])
	plt.xticks([])
	plt.title('Distance to Nearest City')
	plt.savefig(wdfigs+'closestDist',dpi=700)
	
	popClosestDistCities=np.zeros(shape=(iclosestDistCities.shape))
	for ilat in range(len(iclosestDistCities[:,0])):
		for ilon in range(len(iclosestDistCities[0,:])):
			if imageMask2[ilat,ilon]==False:
				popClosestDistCities[ilat,ilon]=cityPop[int(iclosestDistCities[ilat,ilon])]
	np.save(wdvars+'nigeria_popClosestDistCities.npy',popClosestDistCities)

popClosestDistCities=np.ma.masked_array(popClosestDistCities,imageMask2)
closestDistCities=np.ma.masked_array(closestDistCities,imageMask2)
iclosestDistCities=np.ma.masked_array(iclosestDistCities,imageMask2)
plt.clf()
plt.imshow(popClosestDistCities,cmap=cm.viridis_r,vmax=1000000)
plt.colorbar()
plt.yticks([])
plt.xticks([])
plt.title('Population of Nearest City')
plt.savefig(wdfigs+'popclosestDist',dpi=700)


###### Smooth popClosestDistCities ######
try:
	popClosestDistCitiesS=np.load(wdvars+'nigeria_popClosestDistCitiesSmoothed.npy')
	popClosestDistCitiesS=np.ma.masked_array(popClosestDistCitiesS,imageMask2)
except:
	print 'calculating popClosestDistCities Smoothed: try command failed'
	popClosestDistCities1=popClosestDistCities.filled(np.NaN)
	popClosestDistCitiesFilled = replace_nans(popClosestDistCities1, 10, 0.5, 2, method='localmean')
	popClosestDistCitiesS=np.log(popClosestDistCitiesFilled)
	popClosestDistCitiesS=moving_average_2d(popClosestDistCitiesS,np.ones((12, 12)))
	popClosestDistCitiesS=np.ma.masked_array(popClosestDistCitiesS,imageMask2)
	popClosestDistCitiesSmoothed=popClosestDistCitiesS
	popClosestDistCitiesSmoothed[imageMask2==True]=0
	popClosestDistCitiesSmoothed=np.array(popClosestDistCitiesSmoothed)
	np.save(wdvars+'nigeria_popClosestDistCitiesSmoothed',popClosestDistCitiesSmoothed)
popClosestDistCitiesS=scale(popClosestDistCitiesS)+0.01

###### Pop closest city / distance ######
popOverDistCities=np.log(popClosestDistCities/closestDistCities)
popOverDistCities=scale(popOverDistCities)
popClosestDistCitiesS=np.ma.masked_array(popClosestDistCitiesS,imageMask2)

plt.clf()
plt.imshow(popClosestDistCitiesS,cmap=cm.jet)
plt.title('Population of Closest City')
plt.colorbar()
plt.yticks([])
plt.xticks([])
plt.savefig(wdfigs+'popClosestDistCitiesS',dpi=700)
	
plt.clf()
plt.imshow(popOverDistCities,cmap=cm.jet)
plt.title('Population of Closest City / Distance')
plt.yticks([])
plt.xticks([])
plt.colorbar()
plt.savefig(wdfigs+'popOverDistCities',dpi=700)

###########################################
# Night Lights
###########################################
print 'Night Lights'
try:
	stable=np.load(wdvars+'nigeria_stable.npy')
	stableS=np.load(wdvars+'nigeria_stableSmoothed.npy')
	avgVis=np.load(wdvars+'nigeria_avgVis.npy')
	avgVisS=np.load(wdvars+'nigeria_avgVisS.npy')
	avgVisLog=np.load(wdvars+'nigeria_avgVisLog.npy')
	avgVisLogS=np.load(wdvars+'nigeria_avgVisLogS.npy')
except:
	print 'running Night Lights: Try command failed'
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
	
	## grid
	gridl = np.zeros(shape=(len(lonl),len(latl),2))
	for x in range(len(lonl)):
		gridl[x,:,0]=lonl[x]
	for y in range(len(latl)):
		gridl[:,y,1]=latl[y]
	
	## Plot old and new grid
	plt.clf()
	m=4
	s=1152
	plt.plot(np.ma.compressed(gridl[:15,:15,0]),np.ma.compressed(gridl[:15,:15,1]),'*k')
	plt.plot(np.ma.compressed(grid[:m,s-1:s+5,0]),np.ma.compressed(grid[:m,s-1:s+5,1]),'*r')
	plt.yticks([])
	plt.xticks([])
	plt.title('Original Grid (black) and Resized Grid (red)')
	plt.gca().set_aspect('equal', adjustable='box')
	plt.savefig(wdfigs+'old_new_grid.pdf')
	##
	
	stableSmall=stableTif.read_image()
	avgVisSmall=avgVisTif.read_image()
		
	latl=np.radians(latl)
	lonl=np.radians(lonl)
	lutAvgVis=RectSphereBivariateSpline(latl, lonl, avgVisSmall)
	lutStable=RectSphereBivariateSpline(latl, lonl, stableSmall)
	
	newLats=np.radians(latm[::-1])
	newLons=np.radians(lonm)
	newLats,newLons=np.meshgrid(newLats,newLons)
	avgVisRaw=lutAvgVis.ev(newLats.ravel(),newLons.ravel()).reshape((len(lonm),len(latm))).T
	stableRaw=lutStable.ev(newLats.ravel(),newLons.ravel()).reshape((len(lonm),len(latm))).T
	
	
	avgVis=scale(avgVisRaw)
	#try:
	#	distToAvgVis=np.load(wdvars+'nigeria_distToAvgVis.npy')
	#except:
	#	avgVisDist=np.zeros(shape=(avgVis.shape))
	#	avgVis0=np.swapaxes(avgVis,0,1)
	#	avgVis0[avgVis0!=0]=1
	#	avgVisPoints=gridMid[avgVis0[:,:]==1,:]
	#	distToAvgVis,iclosest=findNearest(avgVisPoints,gridMid,imageMask2)
	#	np.save(wdvars+'nigeria_distToAvgVis.npy',distToAvgVis)
	#distToAvgVis=np.ma.masked_array(distToAvgVis,imageMask2)
	
	avgVisS=moving_average_2d(avgVis,np.ones((10, 10)))
	
	avgVisM=np.ma.masked_array(avgVis,imageMask2)
	avgVisLog=np.log(avgVis)
	avgVisLog[np.isinf(avgVisLog)]=np.amin(avgVisLog[~np.isinf(avgVisLog)])
	
	plt.clf()
	plt.imshow(avgVisLog,cmap=cm.jet)
	plt.colorbar()
	plt.title('Avg Visual Lights, log smoothed')
	plt.savefig(wdfigs+'nigeria_avgVisLogS',dpi=700)
	
	avgVisLog=np.clip(scale(avgVisLog),0.94,1)
	avgVisLogS=moving_average_2d(avgVisLog,np.ones((5, 5)))
	
	stable=scale(stableRaw)
	#try:
	#	distToStable=np.load(wdvars+'nigeria_distToStable.npy')
	#except:
	#	stableDist=np.zeros(shape=(stable.shape))
	#	stable0=np.swapaxes(stable,0,1)
	#	stable0[stable0!=0]=1
	#	stablePoints=gridMid[stable0[:,:]==1,:]
	#	distToStable,iclosest=findNearest(stablePoints,gridMid,imageMask2)
	#	np.save(wdvars+'nigeria_distToStable.npy',distToStable)
	#distToStable=np.ma.masked_array(distToStable,imageMask2)
	
	stableS=moving_average_2d(stable,np.ones((10, 10)))
	
	np.save(wdvars+'nigeria_stable.npy',stable)
	np.save(wdvars+'nigeria_stableSmoothed.npy',stableS)
	np.save(wdvars+'nigeria_avgVis.npy',avgVis)
	np.save(wdvars+'nigeria_avgVisS.npy',avgVisS)
	np.save(wdvars+'nigeria_avgVisLog.npy',np.array(avgVisLog))
	np.save(wdvars+'nigeria_avgVisLogS.npy',avgVisLogS)

stable=scale(np.ma.masked_array(stable,imageMask2))
stableS=scale(np.ma.masked_array(stableS,imageMask2))
avgVis=scale(np.ma.masked_array(avgVis,imageMask2))
avgVisLogS=scale(np.ma.masked_array(avgVisLogS,imageMask2))
avgVisS=scale(np.ma.masked_array(avgVisS,imageMask2))

#plt.clf()
#plt.imshow(stable,cmap=cm.jet)
#plt.yticks([])
#plt.xticks([])
#plt.colorbar()
#plt.title('Stable Lights')
#plt.savefig(wdfigs+'nigeria_stable',dpi=700)
#
#plt.clf()
#plt.imshow(stableS,cmap=cm.jet)
#plt.colorbar()
#plt.title('Stable Lights Smoothed')
#plt.savefig(wdfigs+'nigeria_stable_lights',dpi=700)
#
#plt.clf()
#plt.imshow(avgVisLogS,cmap=cm.jet)
#plt.colorbar()
#plt.title('Avg Visual Lights, log smoothed')
#plt.savefig(wdfigs+'nigeria_avgVisLogS',dpi=700)

plt.clf()
plt.imshow(avgVis,cmap=cm.jet)
plt.colorbar()
plt.title('Avg Visual Lights')
plt.savefig(wdfigs+'nigeria_avgVis',dpi=700)

###########################################
# Rivers
###########################################
print 'Rivers'
#### Rivers ####
try:
	riverHits=np.load(wdvars+'nigeria_riverHits.npy')
except:
	print 'calculating riverHits: try command failed'
	rRivers = shapefile.Reader(wddata+'nigeria_landtype/nga_riverspolygon/nga_hyd_riverspolygon_dcw.shp')
	
	records=rRivers.shapeRecords()
	records=MaskableList(records)
	allYear=[rec.record[3] for rec in records]
	waterType=[rec.record[2] for rec in records]
	
	maskPermanent=np.array([rType=='Perennial/Permanent' for rType in allYear])
	maskNonPermanent=np.array([rType=='Non-Perennial/Intermittent/Fluctuating' for rType in allYear])
	maskInlandWater=np.array([rType=='Inland Water' for rType in waterType])
	maskFlood=np.array([rType=='Land Subject to Inundation' for rType in waterType])
	
	permanent=records[maskPermanent]
	nonPermanent=records[maskNonPermanent]
	inland=records[maskInlandWater]
	flood=records[maskFlood]
	
	riverHits=np.zeros(shape=(len(latm),len(lonm),4))
	itype=-1
	for rtype in ['permanent','nonPermanent','inland','flood']:
		itype+=1
		for i in range(len(vars()[rtype])):
			shapePointsLonLat=vars()[rtype][i].shape.points
			shapePointsLonLat=np.array(shapePointsLonLat)
			shapePoints=np.array(shapePointsLonLat)
			shapePoints[:,0]=shapePointsLonLat[:,1]
			shapePoints[:,1]=shapePointsLonLat[:,0]
			riverHits[:,:,itype]+=findGridDataLays(shapePoints,gridMid)
	
	riverHits[riverHits>0]=1
	
	np.save(wdvars+'nigeria_riverHits',riverHits)

riverHits=np.ma.masked_array(riverHits,imageMask3[:,:,:4])

try:
	distToRivers=np.load(wdvars+'nigeria_distToRivers.npy')
except:
	print 'calculating distToRivers: try command failed'
	rivers0=riverHits
	rivers0=1*rivers0
	
	distToRivers=np.zeros(shape=(len(latm),len(lonm),4))
	for iriver in range(4): #'permanent','nonPermanent','inland','flood'
		riversMidPoints=gridMid[rivers0[:,:,iriver]==1,:]
		distToRivers[:,:,iriver],iclosest=findNearest(riversMidPoints,gridMid,imageMask2)
	np.save(wdvars+'nigeria_distToRivers',np.array(distToRivers))
	
distToRivers=np.ma.masked_array(distToRivers,imageMask3[:,:,:len(distToRivers[0,0,:])])

v=-1
for rtype in ['permanent','nonPermanent','inland','flood']:
	v+=1
	plt.clf()
	#plt.imshow(riverHits[:,:,v],cmap=cm.gist_heat_r)
	plt.imshow(distToRivers[:,:,v],vmax=150)
	plt.colorbar()
	plt.title(rtype+' River Density')
	plt.savefig(wdfigs+'nigeria_'+rtype+'_river_density',dpi=700)

###########################################
# Coast Lines
###########################################
print 'Coasts'
try:
	distToCoasts=np.load(wdvars+'nigeria_distToCoasts.npy')
	coastLines=np.load(wdvars+'nigeria_coastLines.npy')
except:
	print 'calculating distToCoasts: try command failed'
	rCoasts0=shapefile.Reader(wddata+'nigeria_landtype/nigeria_coastline/noaa_gshhg/GSHHS_shp/f/GSHHS_f_L1.shp')
	rCoasts1=shapefile.Reader(wddata+'nigeria_landtype/nigeria_coastline/noaa_gshhg/GSHHS_shp/f/GSHHS_f_L2.shp')
	shapes0=rCoasts0.shapes
	shapes1=rCoasts1.shapes
	shapes0=shapes0()
	shapes1=shapes1()
	
	coastLines=np.zeros(shape=(len(latm),len(lonm),2))
	distToCoasts=np.zeros(shape=(len(latm),len(lonm),2))
	for coastType in range(2):
		print '\n Coast',coastType
		for i in range(len(vars()['shapes'+str(coastType)])):
			print round(100*i/float(len(vars()['shapes'+str(coastType)])),2),'%'
			pointsLonLat=np.array(vars()['shapes'+str(coastType)][i].points)
			points=np.array(pointsLonLat)
			points[:,0]=pointsLonLat[:,1]
			points[:,1]=pointsLonLat[:,0]
		
			if (points[:,1]>lonm[0]).any() and (points[:,1]<lonm[-1]).any():
				if (latm[-1]<points[:,1]).any() and (points[:,1]<latm[0]).any():
					# lon mask
					moreThanLon=lonm[0]<points[:,1]
					lessThanLon=points[:,1]<lonm[-1]
					lonMask=moreThanLon==lessThanLon
					# lat mask
					moreThanLat=latm[-1]<points[:,0]
					lessThanLat=points[:,0]<latm[0]
					latMask=moreThanLat==lessThanLat
					# total mask
					nigeriaMask=np.zeros(shape=(len(latMask)),dtype=bool)
					for i in range(len(latMask)):
						if latMask[i]==True and lonMask[i]==True:
							nigeriaMask[i]=True
					if np.sum(nigeriaMask)==0:
						continue
		
					pointsM=points[nigeriaMask]
					coastLines[:,:,coastType]+=findGridDataLays(pointsM,gridMid)
		coastLines[coastLines>1]=1
	
		coastsMidPoints=gridMid[coastLines[:,:,coastType]==1,:]
		distToCoasts[:,:,coastType],iclosest=findNearest(coastsMidPoints,gridMid,imageMask2)
	
	np.save(wdvars+'nigeria_coastLines.npy',coastLines)
	np.save(wdvars+'nigeria_distToCoasts.npy',distToCoasts)

distToCoasts=np.ma.masked_array(distToCoasts,imageMask3[:,:,:2])
coastLines=np.ma.masked_array(coastLines,imageMask3[:,:,:2])

plt.clf()
plt.xticks([])
plt.imshow(np.clip(distToCoasts[:,:,0],0,300))
plt.colorbar()
plt.xticks([])
plt.yticks([])
plt.title('Distance to Coastlines, km')
plt.savefig(wdfigs+'distToCoasts'+str(0),dpi=700)

if makePlots:

	plt.clf()
	map = Basemap(llcrnrlon=lonm[0]-.5,llcrnrlat=np.amin(latm)-.5,urcrnrlon=lonm[-1],urcrnrlat=np.amax(latm), projection='lcc',lon_0=(lonm[-1]+lonm[0])/2,lat_0=(latm[-1]+latm[0])/2,resolution='i')
	map.drawcountries()
	map.readshapefile(wddata+'nigeria_landtype/nigeria_coastline/noaa_gshhg/GSHHS_shp/h/GSHHS_h_L1', name='coastline', drawbounds=True,linewidth=1.5)
	map.readshapefile(wddata+'nigeria_landtype/nigeria_coastline/noaa_gshhg/GSHHS_shp/h/GSHHS_h_L2', name='coastline', drawbounds=True,linewidth=1.5)
	ax = plt.gca()
	plt.title('Nigerian Coastline')
	plt.savefig(wdfigs+'coastline_nigeria_noaa1.pdf')

	coastMask=coastLines==0
	coastLinesM=np.ma.masked_array(coastLines,coastMask)
	plt.clf()
	plt.imshow(imageMask2,cmap=cm.Pastel1)
	plt.imshow(coastLinesM[:,:,0])
	plt.imshow(coastLinesM[:,:,1])
	plt.savefig(wdfigs+'coastLines',dpi=700)

	plt.clf()
	plt.plot(points[latMask,0],points[latMask,1],'*')
	plt.plot(points[lonMask,0],points[lonMask,1],'*')
	plt.plot(points[nigeriaMask,0],points[nigeriaMask,1],'*')
	plt.gca().set_aspect('equal', adjustable='box')
	plt.savefig(wdfigs+'test_coastlines.pdf')


###########################################
# Land Classification
###########################################
print 'Land Classification'
try:
	landCover=np.load(wdvars+'nigeria_landcover.npy')
except:
	ds=gdal.Open(wddata+'nigeria_landtype/MODIS_land_cover/2016_01_01.LC_Type3.tif')
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
		
	landCoversmall=ds.ReadAsArray()
	#landCoversmall=np.ma.masked_array(landCoversmall,landCoversmall==255)
	
	wherex=np.where(landCoversmall==255)[0]
	wherey=np.where(landCoversmall==255)[1]
	for i in range(len(wherex)):
		x=wherex[i]
		y=wherey[i]
		around=landCoversmall[x-2:x+3,y-2:y+3][landCoversmall[x-2:x+3,y-2:y+3]!=255]
		if len(around)==0:
			continue
		landCoversmall[x,y]=mode(around)[0][0]
			
	
	plt.clf()
	plt.imshow(landCoversmall)
	plt.colorbar()
	plt.savefig(wdfigs+'landCoversmall',dpi=700)
	
	landCoverDict={
		0:'Water Bodies',
		1:'Grasslands',
		2:'Shrublands',
		3:'Broadleaf Croplands',
		4:'Savannas',
		5:'Evergreen Broadleaf Forests',
		6:'Deciduous Broadleaf Forests',
		7:'Evergreen Needleleaf Forests',
		8:'Deciduous Needleleaf Forests',
		9:'Non-Vegetated Lands',
		10:'Urban and Built-up Lands'}
	
	### Convert to 5km
	latt=np.radians(latl)
	lont=np.radians(lonl)
	lut=RectSphereBivariateSpline(latt, lont, landCoversmall)
	
	newLats,newLons=np.meshgrid(np.radians(latm[::-1]),np.radians(lonm))
	landCover=lut.ev(newLats.ravel(),newLons.ravel()).reshape((len(lonm),len(latm))).T
	
	###
	
	
	np.save(wdvars+'nigeria_landCover',landCover)

landCover=np.ma.masked_array(landCover,imageMask2)
landCover=np.round(landCover,0)

plt.clf()
plt.imshow(landCover,cmap=cm.gist_ncar,vmin=0)
plt.colorbar()
plt.xticks([])
plt.yticks([])
plt.title('Nigeria Land Cover')
plt.savefig(wdfigs+'landCover',dpi=700)

###########################################
# Trading Indicies
###########################################
print 'Trade'
countries=['Benin','Niger','Chad','Cameroon']
countries=np.array(countries)
countriesDict={'Benin':0,'Niger':1,'Chad':2,'Cameroon':3}

f=open(wddata+'africa_trading_index/gdp_per_capita.csv','r')
countryToGDP={}
gdpPerCapita=np.zeros(shape=len(countries))
i=-2
for line in f:
	i+=1
	line=line.replace('"','')
	tmp=line.split(',')
	if i==-1:
		y2010=np.where(np.array(tmp)=='2010')[0][0]
		continue
	try:
		countryToGDP[tmp[0]]=float(tmp[y2010])
	except:
		continue
	if (tmp[0]==countries).any():
		gdpPerCapita[countriesDict[tmp[0]]]=tmp[y2010]

f=open(wddata+'africa_trading_index/trading_across_borders_allyears.csv','r')
trading=np.zeros(shape=(len(countries),2))
i=-2
for line in f:
	i+=1
	if i==-1:
		continue
	tmp=line.split(',')
	if (tmp[0]==countries).any():
		if tmp[1][2:]=='2010':
			#print tmp[0],tmp[6],tmp[9],float(tmp[6])/float(tmp[9])
			trading[countriesDict[tmp[0]],:]=float(tmp[6]),float(tmp[9])
	
######## Distance to Other Countries ######## 
try:
	borders=np.load(wdvars+'nigeria_borders.npy')
	distToBorders=np.load(wdvars+'nigeria_distToBorders.npy')
except:
	print 'calculating borders: try command failed'
	rCountry=shapefile.Reader(wddata+'africa_trading_index/african_countries_shape/Data/Africa.shp')
	
	records=rCountry.shapeRecords()
	records=MaskableList(records)
	countryID=[rec.record[3] for rec in records]
	
	plt.clf()
	borders=np.zeros(shape=(len(latm),len(lonm),len(countries)))
	k=-1
	for country in countries:
		k+=1
		vars()['mask'+country]=np.array([rType==country for rType in countryID])
		
		vars()['records'+country]=np.array(records[vars()['mask'+country]])
		pointsLonLat=np.array(vars()['records'+country][0].shape.points)
		points=np.array(pointsLonLat)
		points[:,0]=pointsLonLat[:,1]
		points[:,1]=pointsLonLat[:,0]
	
		# lon mask
		moreThanLon=lonm[0]<points[:,1]
		lessThanLon=points[:,1]<lonm[-1]
		lonMask=moreThanLon==lessThanLon
		# lat mask
		moreThanLat=latm[-1]<points[:,0]
		lessThanLat=points[:,0]<latm[0]
		latMask=moreThanLat==lessThanLat
		# total mask
		nigeriaMask=np.zeros(shape=(len(latMask)),dtype=bool)
		for i in range(len(latMask)):
			if latMask[i]==True and lonMask[i]==True:
				nigeriaMask[i]=True
		if np.sum(nigeriaMask)==0:
			continue
		
		pointsM=points[nigeriaMask]
	
		borders[:,:,k]=findGridDataLays(pointsM,gridMid)
		np.save(wdvars+'nigeria_borders.npy',borders)

	distToBorders=np.zeros(shape=(len(latm),len(lonm),len(borders[0,0,:])))
	for iborder in range(len(borders[0,0,:])):
		bordersMidPoints=gridMid[borders[:,:,iborder]==1,:]
		distToBorders[:,:,iborder],iclosest=findNearest(bordersMidPoints,gridMid,imageMask2)
	
		plt.clf()
		plt.imshow(np.clip(distToBorders[:,:,iborder],0,300))
		plt.colorbar()
		plt.title('Distance to '+countries[iborder]+' Border, km (clipped)')
		plt.xticks([])
		plt.yticks([])
		plt.gca().set_aspect('equal', adjustable='box')
		plt.savefig(wdfigs+'distTo'+countries[iborder]+'Border',dpi=700)
	
	np.save(wdvars+'nigeria_distToBorders.npy',distToBorders)

distToBorders=np.ma.masked_array(distToBorders,imageMask3[:,:,:4])
distToBorders=np.clip(distToBorders,0,300)
gdpPerCapitaDistAllCountries=scale(distToBorders*(-gdpPerCapita))
tradeDistAllCountries=np.ma.zeros(shape=(len(latm),len(lonm),4,2))
tradeDistAllCountries[:,:,:,0]=scale(distToBorders*(-trading[:,0]))
tradeDistAllCountries[:,:,:,1]=scale(distToBorders*(-trading[:,1]))

gdpPerCapitaDist=np.sum(gdpPerCapitaDistAllCountries,axis=-1)
tradeDist=np.sum(tradeDistAllCountries,axis=-2)
tradeDist=tradeDist[:,:,1]

plt.clf()
plt.imshow(gdpPerCapitaDist,cmap=cm.jet)
plt.colorbar()
plt.title('GDP of Surrounding Countries by Distance')
plt.gca().set_aspect('equal', adjustable='box')
plt.xticks([])
plt.yticks([])
plt.savefig(wdfigs+'gdpPerCapitaDist',dpi=700)

plt.clf()
plt.imshow(tradeDist,cmap=cm.jet)
plt.colorbar()
plt.title('Trade Ease with Surrounding Countries and Dist')
plt.xticks([])
plt.yticks([])
plt.gca().set_aspect('equal', adjustable='box')
plt.savefig(wdfigs+'tradeDist',dpi=700)
	
bordersInfluence=gdpPerCapitaDist*tradeDist
plt.clf()
plt.imshow(bordersInfluence,cmap=cm.jet)
plt.colorbar()
plt.title('Influence of Borders (GDP and Ease of Trading)')
plt.gca().set_aspect('equal', adjustable='box')
plt.savefig(wdfigs+'bordersInfluence',dpi=700)
	

if makePlots:
	plt.clf()
	masktmp=borders==0
	bordersM=np.ma.masked_array(borders,masktmp)
	plt.imshow(bordersM[:,:,0],cmap=cm.binary,vmin=0,vmax=1)
	plt.imshow(bordersM[:,:,1],cmap=cm.binary,vmin=0,vmax=1)
	plt.imshow(bordersM[:,:,2],cmap=cm.binary,vmin=0,vmax=1)
	plt.imshow(bordersM[:,:,3],cmap=cm.binary,vmin=0,vmax=1)
	plt.colorbar()
	plt.title('Nigerian Borders')
	plt.xticks([])
	plt.yticks([])
	plt.savefig(wdfigs+'nigeria_borders',dpi=700)

	plt.clf()
	map = Basemap(llcrnrlon=lonm[0],llcrnrlat=np.amin(latm),urcrnrlon=lonm[-1],urcrnrlat=np.amax(latm), projection='lcc',lon_0=(lonm[-1]+lonm[0])/2,lat_0=(latm[-1]+latm[0])/2,resolution='i')
	map.readshapefile(wddata+'africa_trading_index/african_countries_shape/Data/Africa', name='countries', drawbounds=True)
	ax = plt.gca()
	plt.title('Country Boundaries')
	plt.savefig(wdfigs+'countries',dpi=700)

	#distToBordersC=distToBordersA
	#distToBordersC[distToBordersC==1000]=0
	#distToBordersBool=1-(distToBordersC==0)
	#distToBordersCount=np.sum(distToBordersBool,axis=-1)
	#distToBordersCount=np.sum(distToBordersC[distToBordersMask],axis=-1)
	plt.clf()
	plt.imshow(gdpPerCapitaDist,cmap=cm.jet_r)
	plt.colorbar()
	plt.title('gdpPerCapita Distance')
	plt.gca().set_aspect('equal', adjustable='box')
	plt.savefig(wdfigs+'gdpPerCapitaDist',dpi=700)
	
	plt.clf()
	plt.imshow(tradingDist,cmap=cm.jet_r)
	plt.colorbar()
	plt.title('Trading and Distance')
	plt.gca().set_aspect('equal', adjustable='box')
	plt.savefig(wdfigs+'gdpPerCapitaDist',dpi=700)
	
###########################################
# Religion
###########################################
#rStates=shapefile.Reader(wddata+'nigeria_landtype/NIR-level_1.shp')
#rRoads = shapefile.Reader(wddata+'nigeria_roads/part/nga_trs_roads_osm_wfpschema_sep2016/nga_trs_roads_osm_wfpschema_sep2016.shp')
#
#rStates=shapefile.Reader(wddata+'boundaries/MapLibraryCD/Data/Africa.shp')
#
#records=rStates.shapeRecords()
#records=MaskableList(records)
#stateName=[rec.record[3] for rec in records1]
#mask1=np.array([Type!='' for Type in stateName])
#records1=records[mask1]
#stateName=[rec.record[3] for rec in records1]
#
#for i in range(len(records1)):
#
#shapePoints=np.array(records1[0].shape.points)
#state,stepsize=findGridDataLays(shapePoints,gridsubsaharan)
#stateFilled=np.zeros(shape=(state.shape))
#
#stateLons=np.where(state==1)[0]
#stateLats=np.where(state==1)[1]
#minLon,minLat=np.amin(stateLons),np.amin(stateLats)
#maxLon,maxLat=np.amax(stateLons),np.amax(stateLats)
#avgLon1=minLon+2
#avgLat2=minLat+2
#whereMinLon1=np.where(stateLons==minLon)[0][len(np.where(stateLons==minLon)[0])/2]
#whereMinLat2=np.where(stateLats==minLat)[0][len(np.where(stateLats==minLat)[0])/2]
#avgLat1=stateLats[whereMinLon1]
#avgLon2=stateLons[whereMinLat2]
#stateFilled[avgLon1,avgLat1]=3
#stateFilled[avgLon2,avgLat2]=3
#stateFilled=stateFilled[minLon:maxLon,minLat:maxLat]
#state=state[minLon:maxLon,minLat:maxLat]
#
#for j in range(2):
#	for ilon in range(len(state[:,0])):
#		if ilon==0:
#			continue
#		for ilat in range(len(state[0,:])):
#			if ilat==0:
#				continue
#			if state[ilon,ilat]==1:
#				continue
#			for k in range(len(state[:,0])):
#				if np.amax(stateFilled[ilon-1:ilon+2,ilat-1:ilat+2])>=1 and np.amax(state[ilon,ilat:ilat+2])==0:
#					stateFilled[ilon,ilat]=1
#	for ilon in range(len(state[:,0]))[::-1]:
#		if ilon==0:
#			continue
#		for ilat in range(len(state[0,:]))[::-1]:
#			if ilat==0:
#				continue
#			for k in range(len(state[:,0])):
#				if np.amax(stateFilled[ilon-1:ilon+2,ilat-1:ilat+2])>=1 and np.amax(state[ilon,ilat:ilat+2])==0:
#					stateFilled[ilon,ilat]=1
#
#for ilon in range(len(state[:,0]))[::-1]:
#	if ilon==0:
#		continue
#	for ilat in range(len(state[0,:]))[::-1]:
#		if ilat==0:
#			continue
#		if ilon==len(state[:,0])-1 or ilat==len(state[0,:])-1:
#			continue
#		if stateFilled[ilon,ilat-1]>=1:
#			stateFilled[ilon,ilat]=1
#		if stateFilled[ilon-1,ilat]>=1:
#			stateFilled[ilon,ilat]=1
#		if np.sum([state[ilon,ilat-1],state[ilon,ilat+1],state[ilon-1,ilat],state[ilon+1,ilat]])>=3:
#			stateFilled[ilon,ilat]=1
#state=np.array(state,dtype=int)
#whereBound=np.where(state==1)
#for k in range(len(whereBound[0])):
#	stateFilled[whereBound[0][k],whereBound[1][k]]=2
#plt.clf()
#plt.imshow(stateFilled)
#plt.savefig(wdfigs+'nigeria_states_filled',dpi=700)
#
#

	#		pointIn=point_in_polygon([ilon,ilat],poly)
	#		if pointIn:
	#			state[ilon,ilat]=1



#	
#	for j in range(len(stateLats)):
#		lon=stateLons[j]
#		lat=stateLats[j]
#		i_other_lats=np.where(stateLats==lat)[0]
#		other_lats=stateLats[i_other_lats]
#		other_lons=stateLons[i_other_lats]
#
#		if len(other_lats)==2 or len(other_lons)==3:
#			state[np.amin(other_lons):np.amax(other_lons),lat]=1
#			continue
#		else:
#			if other_lons
#
#
#		#	other_lons=np.sort(other_lons)
#		#	for i in range(len(other_lats)-1):
#		#		if other_lons[i]==other_lons[i+1] and other_lons[i+1]!=np.amax(other_lons):
#		#				state[other_lons[i]:np.amax(other_lons),lat]=1
#		
#		#if len(other_lats)%2==0:
#		#	other_lons=np.sort(other_lons)
#		#	for k in range(0,len(other_lats),2):
#		#		state[other_lons[k]:other_lons[k+1],lat]=1
#		#	
#		#minLon=np.amax(other_lons)
#		#iminLon=np.where(other_lons==minLon)
#		#lon=stateLons[np.where(stateLats==lat)[0][iminLon][0]]
#		#lat=stateLats[np.where(stateLats==lat)[0][iminLon][0]]
#		#for a in range(1,100):
#		#	if state[lon,lat-a]==1:
#		#		break
#		#	state[lon,lat-a]=1
#		#	exit()
#	plt.clf()
#	plt.imshow(state)
#	plt.savefig(wdfigs+'nigeria_states',dpi=700)
#	exit()
#
#
#
#	exit()
#	plt.plot(shapePoints[:,0],shapePoints[:,1])
#plt.savefig(wdfigs+'nigeria_states.pdf')

#plt.clf()
#map = Basemap(llcrnrlon=-19,llcrnrlat=-35,urcrnrlon=52,urcrnrlat=39, projection='lcc',lon_0=(-19+52)/2,lat_0=(-35+39)/2,resolution='c')
#map.drawcountries()
#map.readshapefile(wddata+'boundaries/MapLibraryCD/Data/Africa', name='admin', drawbounds=True,linewidth=1.5)
#ax = plt.gca()
#plt.title('States')
#plt.savefig(wdfigs+'countires.pdf')

###########################################
# Make Index
###########################################
from sklearn.svm import SVR
import sklearn

############ Multivariate ##############
indicies=[roadDensity[:,:,2],distToRoads[:,:,0],distToRoads[:,:,1],distToRoads[:,:,2],closestDistCities,popOverDistCities,avgVis,distToCoasts[:,:,0],distToCoasts[:,:,1],tradeDist,gdpPerCapitaDist,distToRivers[:,:,1],np.clip(pop5km,0,40000),landCover]
#indicies=[distToRoads[:,:,0],distToRoads[:,:,1],distToRoads[:,:,2],landCover,np.clip(pop5km,0,40000)]
indicies=[np.clip(pop5km,0,40000),distToCoasts[:,:,0]]
xMulti=np.ma.zeros(shape=(len(latm),len(lonm),len(indicies)))
ydata=np.ma.compressed(mal)

for i in range(len(indicies)):
	xMulti[:,:,i]=indicies[i]

multiIndicies=np.array(['roadDensityT','distToRoadsP','distToRoadsS','distToRoadsT','closestDistCities','popClosestDistCitiesS','popOverDistCities','avgVis','distToCoasts','distToInlandCoasts','tradeDist','gdpPerCapitaDist','distToRiversNonPermanent','pop5km','landCover'])

xMultiC=np.zeros(shape=(len(np.ma.compressed(xMulti[:,:,0])),len(xMulti[0,0,:])))
for i in range(len(xMulti[0,0,:])):
	xMulti[:,:,i]=scale(xMulti[:,:,i])
	xMultiC[:,i]=np.ma.compressed(xMulti[:,:,i])

xMultiTrain, xMultiTest, ydataTrain, ydataTest = sklearn.model_selection.train_test_split(xMultiC,ydata,test_size=.1,train_size=.8)

clf=RandomForestRegressor()
clf.fit(xMultiTrain,ydataTrain)

malPred=clf.predict(xMultiTest)
Corr=corr(malPred,ydataTest)
print Corr

x,y=ydataTest*100,malPred*100
slope,b=np.polyfit(x,y,1)
yfit=slope*x+b
error=np.mean((y-x)/x)*100
plt.clf()
plt.plot(x,y,'*b',x,yfit,'g-')
plt.title('Accuracy of Predicted Malnutrition, Corr = '+str(round(Corr,2))+', Avg Error ='+str(round(error,2))+'%')
plt.xlabel('Actual Malnutrition, % population')
plt.ylabel('Predicted Malnutrition, % population')
plt.grid(True)
plt.savefig(wdfigs+'accuracy_predictions.pdf')

#xMulti[:,:,13]=avgVisLog
#xMulti[:,:,13]=avgVisLogS
#xMulti[:,:,14:18]=riverHits
#xMulti[:,:,10]=stable
#xMulti[:,:,11]=stableS
#xMulti[:,:,5]=popClosestDistCities
#xMulti[:,:,5]=popClosestDistCitiesS
# land cover
#xMulti[:,:,22]=landCover
# coasts
#xMulti[:,:,22:24]=distToCoasts
## trading
#xMulti[:,:,24]=tradeDist
#xMulti[:,:,25]=gdpPerCapitaDist
# population

#multiIndicies=['roadDensityP','roadDensityS','roadDensityT','distToRoadsP','distToRoadsS','distToRoadsT','closestDistCities','popClosestDistCities','popClosestDistCitiesS','popOverDistCities','stable','stableS','avgVis','avgVisLogS','riverHitsPermanent','riverHitsNonPermanent','riverHitsInland','riverHitsFlood','distToRiversPermanent','distToRiversNonPermanent','distToRiversHitsInland','distToRiversFlood','distToCoasts','distToInlandCoasts','tradeDist','gdpPerCapitaDist','landCover']
#clf=SVR()
#clf.fit(xMultiC,ydata)
#coef=clf.coef_
#coef=clf.feature_importances_
#sortedCoefIndex=np.argsort(abs(coef))
#sortedCoefIndex=sortedCoefIndex[::-1]
#print multiIndicies[sortedCoefIndex]
#print coef[sortedCoefIndex]

#malPred=xMulti[:,:,:]*coef
#malPred=np.sum(malPred,axis=-1)
malPred=np.zeros(shape=(landCover.shape))
malPred=np.ma.masked_array(malPred,imageMask2)
x1=0
for i in range(len(malPred[:,0])):
	x2=len(malPred[i][np.ma.getmask(malPred[i])==False])+x1
	malPred[i][np.ma.getmask(malPred[i])==False]=clf.predict(xMultiC[x1:x2])
	x1=x2
#malPred=clf.predict(xMultiC)

plt.clf()
plt.imshow(malPred,cmap=cm.jet,vmax=0.3)
plt.colorbar()
plt.yticks([])
plt.xticks([])
plt.title('Predicted Malnutrition, Percent of Population')
plt.savefig(wdfigs+'malPred',dpi=700)
exit()

lutPovPred=RectSphereBivariateSpline(oldLat, oldLon, povPred)
povPredC=lutPovPred.ev(newLats.ravel(),newLons.ravel()).reshape((widthC,heightC)).T
povPredC=np.ma.masked_array(povPredC,maskCoarse)

lutPov125=RectSphereBivariateSpline(oldLat, oldLon, np.array(pov125))
pov125C=lutPov125.ev(newLats.ravel(),newLons.ravel()).reshape((widthC,heightC)).T
pov125C=np.ma.masked_array(pov125C,maskCoarse)
pov125C=np.ma.masked_array(pov125C,np.isnan(pov125C))

plt.clf()
plt.imshow(povPredC,cmap=cm.jet_r)
plt.colorbar()
plt.yticks([])
plt.xticks([])
plt.title('Predicted Poverty, 5km')
plt.savefig(wdfigs+'predPov',dpi=700)

plt.clf()
plt.imshow(pov125C,cmap=cm.jet_r)
plt.colorbar()
plt.yticks([])
plt.xticks([])
plt.title('Percent Living Under 1.25/day 2010, 5km')
plt.savefig(wdfigs+'nigeria_pov125C',dpi=700)

