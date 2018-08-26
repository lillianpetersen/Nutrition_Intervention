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
travelbi=np.array(travel)
travelbi[travel!=0]=1
travelbi=1-travelbi

plt.clf()
plt.imshow(travelbi,cmap=cm.gist_ncar_r)
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

roadDensity=np.swapaxes(roadDensity[:,:,:],0,1)
roadDensity=np.ma.masked_array(roadDensity,imageMask3[:,:,:3])

if makePlots:
	plt.clf()
	map = Basemap(llcrnrlon=lonM[0],llcrnrlat=np.amin(latM),urcrnrlon=lonM[-1],urcrnrlat=np.amax(latM), projection='lcc',lon_0=(lonM[-1]+lonM[0])/2,lat_0=(latM[-1]+latM[0])/2,resolution='i')
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
	
	distToRoads=1-roads0 # raod = 0
	distToRoads[distToRoads==1]=9999.
	
	distToRoads=np.zeros(shape=(len(latM),len(lonM),3))
	for iroad in range(3):
		roadsMidPoints=gridMid[roads0[:,:,iroad]==1,:]
		print 'finding distance to roads',iroad
		distToRoads[:,:,iroad],iclosestRoad=findNearest(roadsMidPoints,gridMid,imageMask1)
	np.save(wdvars+'nigeria_distToRoads',np.array(distToRoads))
	distToRoads=np.ma.masked_array(distToRoads,imageMask3[:,:,:3])
		
if makePlots:
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
	cityLonLat=np.load(wdvars+'NigeriaCityLonLat.npy')
except:
	print 'calculating cityPop: try command failed'
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
	np.save(wdvars+'NigeriaCityPop',cityPop)
	np.save(wdvars+'NigeriaCityLonLat',cityLonLat)

###### Dist To Cities and Pop of Closest City ######
try:
	closestDistCities=np.load(wdvars+'closestDist_to_nigerian_cities.npy')
	iclosestDistCities=np.load(wdvars+'iclosestDist_to_nigerian_cities.npy')
	popClosestDistCities=np.load(wdvars+'nigeria_popClosestDistCities.npy')
except:
	print 'calculating closestDistCities: try command failed'
	closestDistCities,iclosestDistCities=findNearest(cityLonLat,gridMid,imageMask1)
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

###### Smooth popClosestDistCities ######
try:
	popClosestDistCitiesS=np.load(wdvars+'nigeria_popClosestDistCitiesSmoothed.npy')
	popClosestDistCitiesS=np.ma.masked_array(popClosestDistCitiesS,imageMask2)
except:
	print 'calculating popClosestDistCities Smoothed: try command failed'
	popClosestDistCities1=popClosestDistCitiesCities.filled(np.NaN)
	popClosestDistCitiesFilled = replace_nans(popClosestDistCities1, 10, 0.5, 2, method='localmean')
	popClosestDistCitiesS=np.log(popClosestDistCitiesFilled)
	popClosestDistCitiesS=moving_average_2d(popClosestDistCitiesS,np.ones((40, 40)))
	popClosestDistCitiesS=np.ma.masked_array(popClosestDistCitiesS,imageMask2)
	popClosestDistCitiesSmoothed=popClosestDistCitiesS
	popClosestDistCitiesSmoothed[imageMask2==True]=0
	popClosestDistCitiesSmoothed=np.array(popClosestDistCitiesSmoothed)
	np.save(wdvars+'nigeria_popClosestDistCitiesSmoothed',popClosestDistCitiesSmoothed)
popClosestDistCitiesS=scale(popClosestDistCitiesS)+0.01

###### Pop closest city / distance ######
popOverDistCities=np.log(popClosestDistCities/closestDistCities)
popOverDistCities=scale(popOverDistCities)
if makePlots:
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
	
	newLats=np.radians(latM[::-1])
	newLons=np.radians(lonM)
	newLats,newLons=np.meshgrid(newLats,newLons)
	avgVisRaw=lutAvgVis.ev(newLats.ravel(),newLons.ravel()).reshape((len(lonM),len(latM))).T
	stableRaw=lutStable.ev(newLats.ravel(),newLons.ravel()).reshape((len(lonM),len(latM))).T
	
	
	avgVis=scale(avgVisRaw)
	#try:
	#	distToAvgVis=np.load(wdvars+'nigeria_distToAvgVis.npy')
	#except:
	#	avgVisDist=np.zeros(shape=(avgVis.shape))
	#	avgVis0=np.swapaxes(avgVis,0,1)
	#	avgVis0[avgVis0!=0]=1
	#	avgVisPoints=gridMid[avgVis0[:,:]==1,:]
	#	distToAvgVis,iclosest=findNearest(avgVisPoints,gridMid,imageMask1)
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
	#	distToStable,iclosest=findNearest(stablePoints,gridMid,imageMask1)
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

if makePlots:
	plt.clf()
	plt.imshow(stable,cmap=cm.jet)
	plt.yticks([])
	plt.xticks([])
	plt.colorbar()
	plt.title('Stable Lights')
	plt.savefig(wdfigs+'nigeria_stable',dpi=700)
	
	plt.clf()
	plt.imshow(stableS,cmap=cm.jet)
	plt.colorbar()
	plt.title('Stable Lights Smoothed')
	plt.savefig(wdfigs+'nigeria_stable_lights',dpi=700)
	
	plt.clf()
	plt.imshow(avgVis,cmap=cm.jet)
	plt.colorbar()
	plt.title('Avg Visual Lights')
	plt.savefig(wdfigs+'nigeria_avgVis',dpi=700)
	
	plt.clf()
	plt.imshow(avgVisLogS,cmap=cm.jet)
	plt.colorbar()
	plt.title('Avg Visual Lights, log smoothed')
	plt.savefig(wdfigs+'nigeria_avgVisLogS',dpi=700)

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
	
	riverHits=np.zeros(shape=(len(lonM),len(latM),4))
	itype=-1
	for rtype in ['permanent','nonPermanent','inland','flood']:
		itype+=1
		for i in range(len(vars()[rtype])):
			shapePoints=vars()[rtype][i].shape.points
			riverHits[:,:,itype]+=findGridDataLays(shapePoints,grid)
	
	riverHits[riverHits>0]=1
	
	riverHits=np.swapaxes(riverHits[:,:,:],0,1)
	np.save(wdvars+'nigeria_riverHits',riverHits)

riverHits=np.ma.masked_array(riverHits,imageMask3[:,:,:4])

try:
	distToRivers=np.load(wdvars+'nigeria_distToRivers.npy')
except:
	print 'calculating distToRivers: try command failed'
	rivers0=riverHits
	rivers0=1*rivers0
	rivers0=np.swapaxes(rivers0,0,1)
	
	#distToRivers=np.zeros(shape=(len(latM),len(lonM),4))
	for iriver in range(4): #'permanent','nonPermanent','inland','flood'
		riversMidPoints=gridMid[rivers0[:,:,iriver]==1,:]
		distToRivers[:,:,iriver],iclosest=findNearest(riversMidPoints,gridMid,imageMask1)
	np.save(wdvars+'nigeria_distToRivers',np.array(distToRivers))
	
distToRivers=np.ma.masked_array(distToRivers,imageMask3[:,:,:len(distToRivers[0,0,:])])

if makePlots:
	v=-1
	for rtype in ['permanent','nonPermanent','inland','flood']:
		v+=1
		plt.clf()
		#plt.imshow(riverHits[:,:,v],cmap=cm.gist_heat_r)
		plt.imshow(riverHits[:,:,v],cmap=cm.nipy_spectral_r)
		plt.colorbar()
		plt.title(rtype+' River Density')
		plt.savefig(wdfigs+'nigeria_'+rtype+'_river_density',dpi=700)

plt.clf()
plt.imshow(np.clip(distToRivers[:,:,1],0,60))
plt.title('Dist To Floodlands')
plt.xticks([])
plt.yticks([])
plt.colorbar()
plt.savefig(wdfigs+'dist_to_nonPermanent_Rivers',dpi=700)

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
	
	coastLines=np.zeros(shape=(len(lonM),len(latM),2))
	distToCoasts=np.zeros(shape=(len(latM),len(lonM),2))
	for coastType in range(2):
		for i in range(len(vars()['shapes'+str(coastType)])):
			points=np.array(vars()['shapes'+str(coastType)][i].points)
		
			if (points[:,0]>lonM[0]).any() and (points[:,0]<lonM[-1]).any():
				if (latM[-1]<points[:,1]).any() and (points[:,1]<latM[0]).any():
					# lon mask
					moreThanLon=lonM[0]<points[:,0]
					lessThanLon=points[:,0]<lonM[-1]
					lonMask=moreThanLon==lessThanLon
					# lat mask
					moreThanLat=latM[-1]<points[:,1]
					lessThanLat=points[:,1]<latM[0]
					latMask=moreThanLat==lessThanLat
					# total mask
					nigeriaMask=np.zeros(shape=(len(latMask)),dtype=bool)
					for i in range(len(latMask)):
						if latMask[i]==True and lonMask[i]==True:
							nigeriaMask[i]=True
					if np.sum(nigeriaMask)==0:
						continue
		
					pointsM=points[nigeriaMask]
					coastLines[:,:,coastType]=findGridDataLays(pointsM,grid)
	
		coastsMidPoints=gridMid[coastLines[:,:,coastType]==1,:]
		distToCoasts[:,:,coastType],iclosest=findNearest(coastsMidPoints,gridMid,imageMask1)
	
	coastLines=np.swapaxes(coastLines,0,1)
	np.save(wdvars+'nigeria_coastLines.npy',coastLines)
	np.save(wdvars+'nigeria_distToCoasts.npy',distToCoasts)

distToCoasts=np.ma.masked_array(distToCoasts,imageMask3[:,:,:2])
coastLines=np.ma.masked_array(coastLines,imageMask3[:,:,:2])

if makePlots:
	plt.clf()
	plt.xticks([])
	plt.imshow(np.clip(distToCoasts[:,:,0],0,300))
	plt.colorbar()
	plt.xticks([])
	plt.yticks([])
	plt.title('Distance to Coastlines, km')
	plt.savefig(wdfigs+'distToCoasts'+str(0),dpi=700)

	plt.clf()
	map = Basemap(llcrnrlon=lonM[0]-.5,llcrnrlat=np.amin(latM)-.5,urcrnrlon=lonM[-1],urcrnrlat=np.amax(latM), projection='lcc',lon_0=(lonM[-1]+lonM[0])/2,lat_0=(latM[-1]+latM[0])/2,resolution='i')
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
	landTypeTif=TIFF.open(wddata+'nigeria_landtype/nigeriaLandClassification.tif',mode='r')
	
	ds=gdal.Open(wddata+'nigeria_landtype/nigeriaLandClassification.tif')
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
	latl=np.radians(latl)
	lonl=np.radians(lonl)
	
	landCoverFull=landTypeTif.read_image()
	landCoverDict={
		11:'Post-flooding or irrigated croplands',
		14:'Rainfed croplands',
		20:'Mosaic cropland (50-70%) / vegetation (grassland shrubland forest) (20-50%)',
		30:'Mosaic vegetation (grassland shrubland forest) (50-70%) / cropland (20-50%)',
		40:'Closed to open (>15%) broadleaved evergreen and/or semi-deciduous forest (>5m)',
		50:'Closed (>40%) broadleaved deciduous forest (>5m)',
		60:'Open (15-40%) broadleaved deciduous forest (>5m)',
		70:'Closed (>40%) needleleaved evergreen forest (>5m)',
		90:'Open (15-40%) needleleaved deciduous or evergreen forest (>5m)',
		100:'Closed to open (>15%) mixed broadleaved and needleleaved forest (>5m)',
		110:'Mosaic forest-shrubland (50-70%) / grassland (20-50%)',
		120:'Mosaic grassland (50-70%) / forest-shrubland (20-50%)',
		130:'Closed to open (>15%) shrubland (<5m)',
		140:'Closed to open (>15%) grassland',
		150:'Sparse (>15%) vegetation (woody vegetation shrubs grassland)',
		160:'Closed (>40%) broadleaved forest regularly flooded - Fresh water',
		170:'Closed (>40%) broadleaved semi-deciduous and/or evergreen forest regularly flooded - saline water',
		180:'Closed to open (>15%) vegetation (grassland shrubland woody vegetation) on regularly flooded or waterlogged soil - fresh brackish or saline water',
		190:'Artificial surfaces and associated areas (urban areas >50%) GLOBCOVER 2009',
		200:'Bare areas',
		210:'Water bodies',
		220:'Permanent snow and ice',
		230:'Unclassified'}
	
	landCoverDictNew={
		11:'Post-flooding or irrigated croplands',
		14:'Rainfed croplands',
		20:'Mosaic cropland (50-70%) / vegetation (grassland shrubland forest) (20-50%)',
		30:'Mosaic vegetation (grassland shrubland forest) (50-70%) / cropland (20-50%)',
		40:'Closed to open (>15%) broadleaved evergreen and/or semi-deciduous forest (>5m)',
		50:'Closed (>40%) broadleaved deciduous forest (>5m)',
		60:'Open (15-40%) broadleaved deciduous forest (>5m)',
		70:'Closed (>40%) needleleaved evergreen forest (>5m)',
		90:'Open (15-40%) needleleaved deciduous or evergreen forest (>5m)',
		100:'Closed to open (>15%) mixed broadleaved and needleleaved forest (>5m)',
		110:'Mosaic forest-shrubland (50-70%) / grassland (20-50%)',
		120:'Mosaic grassland (50-70%) / forest-shrubland (20-50%)',
		130:'Closed to open (>15%) shrubland (<5m)',
		140:'Closed to open (>15%) grassland',
		150:'Sparse (>15%) vegetation (woody vegetation shrubs grassland)',
		160:'Closed (>40%) broadleaved forest regularly flooded - Fresh water',
		170:'Closed (>40%) broadleaved semi-deciduous and/or evergreen forest regularly flooded - saline water',
		180:'Closed to open (>15%) vegetation (grassland shrubland woody vegetation) on regularly flooded or waterlogged soil - fresh brackish or saline water',
		210:'Water bodies',
		190:'Artificial surfaces and associated areas (urban areas >50%) GLOBCOVER 2009',
		200:'Bare areas',
		220:'Permanent snow and ice',
		230:'Unclassified'}
	
	convertLandCover=str(landCoverDictNew).replace('{','').replace('}','').split(',')
	#convertLandCover=str(landCoverConversion).replace('{','').replace('}','').split(',')
	landCoverSorted=np.array(landCoverFull)
	for i in range(len(convertLandCover)):
		landCoverSorted[landCoverSorted==int(convertLandCover[i].split(':')[0])]=i
	
	plt.clf()
	plt.imshow(landCoverSorted,cmap=cm.jet_r)
	plt.colorbar()
	plt.savefig(wdfigs+'landCover_fullres_sorted',dpi=700)
	
	### Convert to 1km
	lutland=RectSphereBivariateSpline(latl, lonl, landCoverSorted)
	
	newLats=np.radians(latM[::-1])
	newLons=np.radians(lonM)
	newLats,newLons=np.meshgrid(newLats,newLons)
	landCover=lutland.ev(newLats.ravel(),newLons.ravel()).reshape((len(lonM),len(latM))).T
	
	plt.clf()
	plt.imshow(landCover,cmap=cm.gist_ncar,vmin=0)
	plt.colorbar()
	plt.xticks([])
	plt.yticks([])
	plt.title('Nigeria Land Cover')
	plt.savefig(wdfigs+'landCover',dpi=700)
	
	np.save(wdvars+'nigeria_landCover',landCover)

landCover=np.ma.masked_array(landCover,imageMask2)

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
	borders=np.zeros(shape=(len(lonM),len(latM),len(countries)))
	k=-1
	for country in countries:
		k+=1
		vars()['mask'+country]=np.array([rType==country for rType in countryID])
		
		vars()['records'+country]=np.array(records[vars()['mask'+country]])
		points=np.array(vars()['records'+country][0].shape.points)
	
		# lon mask
		moreThanLon=lonM[0]<points[:,0]
		lessThanLon=points[:,0]<lonM[-1]
		lonMask=moreThanLon==lessThanLon
		# lat mask
		moreThanLat=latM[-1]<points[:,1]
		lessThanLat=points[:,1]<latM[0]
		latMask=moreThanLat==lessThanLat
		# total mask
		nigeriaMask=np.zeros(shape=(len(latMask)),dtype=bool)
		for i in range(len(latMask)):
			if latMask[i]==True and lonMask[i]==True:
				nigeriaMask[i]=True
		if np.sum(nigeriaMask)==0:
			continue
		
		pointsM=points[nigeriaMask]
	
		borders[:,:,k]=findGridDataLays(pointsM,grid)
		np.save(wdvars+'nigeria_borders.npy',borders)

	distToBorders=np.zeros(shape=(len(latM),len(lonM),len(borders[0,0,:])))
	for iborder in range(len(borders[0,0,:])):
		#bordersMidPoints=gridMid[borders[:,:,iborder]==1,:]
		#distToBorders[:,:,iborder],iclosest=findNearest(bordersMidPoints,gridMid,imageMask1)
	
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
tradeDistAllCountries=np.ma.zeros(shape=(len(latM),len(lonM),4,2))
tradeDistAllCountries[:,:,:,0]=scale(distToBorders*(-trading[:,0]))
tradeDistAllCountries[:,:,:,1]=scale(distToBorders*(-trading[:,1]))

gdpPerCapitaDist=np.sum(gdpPerCapitaDistAllCountries,axis=-1)
tradeDist=np.sum(tradeDistAllCountries,axis=-2)
tradeDist=tradeDist[:,:,1]

borders=np.swapaxes(borders,0,1)


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
	map = Basemap(llcrnrlon=lonM[0],llcrnrlat=np.amin(latM),urcrnrlon=lonM[-1],urcrnrlat=np.amax(latM), projection='lcc',lon_0=(lonM[-1]+lonM[0])/2,lat_0=(latM[-1]+latM[0])/2,resolution='i')
	map.readshapefile(wddata+'africa_trading_index/african_countries_shape/Data/Africa', name='countries', drawbounds=True)
	ax = plt.gca()
	plt.title('Country Boundaries')
	plt.savefig(wdfigs+'countries',dpi=700)

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
rStates=shapefile.Reader(wddata+'nigeria_landtype/NIR-level_1.shp')
#rRoads = shapefile.Reader(wddata+'nigeria_roads/part/nga_trs_roads_osm_wfpschema_sep2016/nga_trs_roads_osm_wfpschema_sep2016.shp')

records=rStates.shapeRecords()
records=MaskableList(records)
stateName=[rec.record[3] for rec in records]
mask1=np.array([Type!='' for Type in stateName])
records1=records[mask1]

for i in range(len(records1)):
	shapePoints=np.array(records1[i].shape.points)
	state=findGridDataLays(shapePoints,grid)
	stateFilled=np.zeros(shape=(state.shape))

	stateLons=np.where(state==1)[0]
	stateLats=np.where(state==1)[1]
	minLon,minLat=np.amin(stateLons),np.amin(stateLats)
	maxLon,maxLat=np.amax(stateLons),np.amax(stateLats)
	avgLon1=minLon+2
	avgLat2=minLat+2
	whereMinLon1=np.where(stateLons==minLon)[0][len(np.where(stateLons==minLon)[0])/2]
	whereMinLat2=np.where(stateLats==minLat)[0][len(np.where(stateLats==minLat)[0])/2]
	avgLat1=stateLats[whereMinLon1]
	avgLon2=stateLons[whereMinLat2]
	stateFilled[avgLon1,avgLat1]=3
	stateFilled[avgLon2,avgLat2]=3
	stateFilled=stateFilled[minLon:maxLon,minLat:maxLat]
	state=state[minLon:maxLon,minLat:maxLat]

	for j in range(2):
		for ilon in range(len(state[:,0])):
			if ilon==0:
				continue
			for ilat in range(len(state[0,:])):
				if ilat==0:
					continue
				if state[ilon,ilat]==1:
					continue
				for k in range(len(state[:,0])):
					if np.amax(stateFilled[ilon-1:ilon+2,ilat-1:ilat+2])>=1 and np.amax(state[ilon,ilat:ilat+2])==0:
						stateFilled[ilon,ilat]=1
		for ilon in range(len(state[:,0]))[::-1]:
			if ilon==0:
				continue
			for ilat in range(len(state[0,:]))[::-1]:
				if ilat==0:
					continue
				for k in range(len(state[:,0])):
					if np.amax(stateFilled[ilon-1:ilon+2,ilat-1:ilat+2])>=1 and np.amax(state[ilon,ilat:ilat+2])==0:
						stateFilled[ilon,ilat]=1

	for ilon in range(len(state[:,0]))[::-1]:
		if ilon==0:
			continue
		for ilat in range(len(state[0,:]))[::-1]:
			if ilat==0:
				continue
			if ilon==len(state[:,0])-1 or ilat==len(state[0,:])-1:
				continue
			if stateFilled[ilon,ilat-1]>=1:
				stateFilled[ilon,ilat]=1
			if stateFilled[ilon-1,ilat]>=1:
				stateFilled[ilon,ilat]=1
			if np.sum([state[ilon,ilat-1],state[ilon,ilat+1],state[ilon-1,ilat],state[ilon+1,ilat]])>=3:
				stateFilled[ilon,ilat]=1
	state=np.array(state,dtype=int)
	whereBound=np.where(state==1)
	for k in range(len(whereBound[0])):
		stateFilled[whereBound[0][k],whereBound[1][k]]=2
	plt.clf()
	plt.imshow(stateFilled)
	plt.savefig(wdfigs+'nigeria_states_filled',dpi=700)
	exit()



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

plt.clf()
map = Basemap(llcrnrlon=lonM[0]-.5,llcrnrlat=np.amin(latM)-.5,urcrnrlon=lonM[-1],urcrnrlat=np.amax(latM), projection='lcc',lon_0=(lonM[-1]+lonM[0])/2,lat_0=(latM[-1]+latM[0])/2,resolution='i')
map.drawcountries()
map.readshapefile(wddata+'nigeria_landtype/NIR-level_1', name='admin', drawbounds=True,linewidth=1.5)
ax = plt.gca()
plt.title('Nigerian States')
plt.savefig(wdfigs+'nigeria_states.pdf')

'''
###########################################
# Make Index
###########################################
############ Multivariate ##############
xMulti=np.ma.zeros(shape=(len(latM),len(lonM),23))
ydata=np.ma.compressed(pov125)
# roads
xMulti[:,:,0:3]=roadDensity
xMulti[:,:,3:6]=distToRoads
# cities
xMulti[:,:,6]=closestDistCities
xMulti[:,:,7]=popClosestDistCities
xMulti[:,:,8]=popClosestDistCitiesS
xMulti[:,:,9]=popOverDistCities
# lights
xMulti[:,:,10]=stable
xMulti[:,:,11]=stableS
xMulti[:,:,12]=avgVis
#xMulti[:,:,13]=avgVisLog
xMulti[:,:,13]=avgVisLogS
# coasts
xMulti[:,:,14:16]=distToCoasts
# trading
xMulti[:,:,16]=tradeDist
xMulti[:,:,17]=gdpPerCapitaDist
# rivers
xMulti[:,:,14:18]=riverHits
xMulti[:,:,18:22]=distToRivers
# land cover
xMulti[:,:,22]=landCover
# coasts
#xMulti[:,:,22:24]=distToCoasts
## trading
#xMulti[:,:,24]=tradeDist
#xMulti[:,:,25]=gdpPerCapitaDist

multiIndicies=['roadDensityP','roadDensityS','roadDensityT','distToRoadsP','distToRoadsS','distToRoadsT','closestDistCities','popClosestDistCities','popClosestDistCitiesS','popOverDistCities','stable','stableS','avgVis','avgVisLogS','riverHitsPermanent','riverHitsNonPermanent','riverHitsInland','riverHitsFlood','distToRiversPermanent','distToRiversNonPermanent','distToRiversHitsInland','distToRiversFlood','distToCoasts','distToInlandCoasts','tradeDist','gdpPerCapitaDist','landCover']
multiIndicies=np.array(multiIndicies)

xMultiC=np.zeros(shape=(len(np.ma.compressed(xMulti[:,:,0])),len(xMulti[0,0,:])))
for i in range(len(xMulti[0,0,:])):
	xMulti[:,:,i]=scale(xMulti[:,:,i])
	xMultiC[:,i]=np.ma.compressed(xMulti[:,:,i])

clf=linear_model.LinearRegression()
clf.fit(xMultiC,ydata)
Corr=np.sqrt(clf.score(xMultiC,ydata))
coef=clf.coef_
sortedCoefIndex=np.argsort(abs(coef))
sortedCoefIndex=sortedCoefIndex[::-1]
print multiIndicies[sortedCoefIndex]
print coef[sortedCoefIndex]

povPred=xMulti[:,:,:]*coef
povPred=np.sum(povPred,axis=-1)
povPred=scale(povPred)

plt.clf()
plt.imshow(povPred*100,cmap=cm.jet_r)
plt.colorbar()
plt.yticks([])
plt.xticks([])
plt.title('Predicted Poverty, % living under 1.25 US$/day')
plt.savefig(wdfigs+'povPred',dpi=700)
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

exit()

#np.save(wdvars+'Illinois/xMulti',xMulti)
#np.save(wdvars+'Illinois/ydataMulti',ydata)
#### For all of Nigeria ####
#distToPRoadsS=(distToRoads[:,:,0]-np.amin(distToRoads[:,:,0]))/(np.amax(distToRoads[:,:,0])-np.amin(distToRoads[:,:,0]))
#distToSRoadsS=(distToRoads[:,:,1]-np.amin(distToRoads[:,:,1]))/(np.amax(distToRoads[:,:,1])-np.amin(distToRoads[:,:,1]))
# Primary
distToRoadsClipped=np.clip(distToRoads[:,:,0],0,15)
distToPRoadsS=(distToRoadsClipped[:,:]-np.amin(distToRoadsClipped[:,:]))/(np.amax(distToRoadsClipped[:,:])-np.amin(distToRoadsClipped[:,:]))
# Secondary
distToRoadsClipped=np.clip(distToRoads[:,:,1],0,10)
distToSRoadsS=(distToRoadsClipped[:,:]-np.amin(distToRoadsClipped[:,:]))/(np.amax(distToRoadsClipped[:,:])-np.amin(distToRoadsClipped[:,:]))
# Teritiary
distToRoadsClipped=np.clip(distToRoads[:,:,2],0,3)
distToTRoadsS=(distToRoadsClipped[:,:]-np.amin(distToRoadsClipped[:,:]))/(np.amax(distToRoadsClipped[:,:])-np.amin(distToRoadsClipped[:,:]))

TroadDensityS=(roadDensity[:,:,2]-np.amin(roadDensity[:,:,2]))/(np.amax(roadDensity[:,:,2])-np.amin(roadDensity[:,:,2]))
SroadDensityS=(roadDensity[:,:,1]-np.amin(roadDensity[:,:,1]))/(np.amax(roadDensity[:,:,1])-np.amin(roadDensity[:,:,1]))
#roadDensityS=popclosestDistCities/roadDensityS
closestDistS=(closestDistCities-0)/(np.amax(closestDistCities)-0)

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

indext1=np.log(popClosestDistCitiesS/closestDistS)
indext1=(indext1-np.amin(indext1))/(np.amax(indext1)-np.amin(indext1))
indext1=1-indext1
#indext1=closestDistS
plt.clf()
plt.imshow(indext1,cmap=cm.hot_r)
plt.colorbar()
plt.title('popCity/DistToCity')
plt.savefig(wdfigs+'indext1',dpi=700)

#index1=distToPRoadsS-roadDensityS+closestDistS-popClosestDistCitiesS
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
oldLat,oldLon=np.radians(latM[::-1]),np.radians(lonM)
lutIndex1=RectSphereBivariateSpline(oldLat, oldLon, index1)
lutPov=RectSphereBivariateSpline(oldLat, oldLon, pov125)
lutMask=RectSphereBivariateSpline(oldLat, oldLon, imageMask2)

newLats,newLons=np.radians(newLats),np.radians(newLons)
newLats,newLons=np.meshgrid(newLats,newLons)

urbanIndex=lutIndex1.ev(newLats.ravel(),newLons.ravel()).reshape((width,height)).T
pov125Coarse=lutPov.ev(newLats.ravel(),newLons.ravel()).reshape((width,height)).T
maskCoarse=lutMask.ev(newLats.ravel(),newLons.ravel()).reshape((width,height)).T
maskCoarse=np.array(np.round(maskCoarse,0),dtype=bool)

urbanIndex=np.ma.masked_array(urbanIndex,maskCoarse)
urbanIndex=(urbanIndex-np.amin(urbanIndex))/(np.amax(urbanIndex)-np.amin(urbanIndex))
pov125Coarse=np.ma.masked_array(pov125Coarse,maskCoarse)

plt.clf()
plt.imshow(urbanIndex,cmap=cm.jet_r)
plt.colorbar()
plt.title('Urban/Rural Index Based on Roads and Population')
plt.savefig(wdfigs+'urbanIndex',dpi=700)

plt.clf()
plt.imshow(pov125Coarse,cmap=cm.jet_r)
plt.colorbar()
plt.title('Percent Under $1.25/Day')
plt.savefig(wdfigs+'pov125Coarse',dpi=700)
print corr(np.ma.compressed(urbanIndex),np.ma.compressed(pov125Coarse))

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
plt.xlabel('distToPRoads+distToSRoads-(TroadDensity+SroadDensity)+log(popclosestDistCities/closestDistCities)')
plt.savefig(wdfigs+'index1_and_pov125_nigeria',dpi=700)
exit()


plt.clf()
plt.imshow(index,cmap=cm.hot_r,norm=colors.LogNorm())
plt.colorbar()
plt.title('Rural Index by Cities')
plt.savefig(wdfigs+'index',dpi=700)
distToAvgVis=np.ma.masked_array(distToAvgVis,imageMask2)


plt.clf()
plt.imshow(iclosestDistCities)
plt.colorbar()
plt.title('Nearest City')
plt.savefig(wdfigs+'iclosestDist',dpi=700)

plt.clf()
plt.imshow(popclosestDistCities)
plt.colorbar()
plt.title('Nearest City Population')
plt.savefig(wdfigs+'popClosestDistCities',dpi=700)

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
'''
