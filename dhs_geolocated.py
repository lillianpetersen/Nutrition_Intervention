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
# Read Poverty
##############################################
country='nigeria'

plt.clf()
map = Basemap(llcrnrlon=lonM[0],llcrnrlat=np.amin(latM),urcrnrlon=lonM[-1],urcrnrlat=np.amax(latM), projection='lcc',lon_0=(lonM[-1]+lonM[0])/2,lat_0=(latM[-1]+latM[0])/2,resolution='i')
map.drawcoastlines(linewidth=1.5)
map.drawcountries()
map.readshapefile(wddata+'dhs_nigeria_geo/NGGE71FL/NGGE71FL', name='pov', drawbounds=True)
ax = plt.gca()
plt.title('Nigerian DHS Clusters')
plt.savefig(wdfigs+'clusters',dpi=700)

rDHS = shapefile.Reader(wddata+'dhs_nigeria_geo/NGGE71FL/NGGE71FL.shp')

records=rDHS.shapeRecords()
points=[rec.shape.points for rec in records]
points=np.array(points)
points=points[:,0,:]
points[np.round(points,4)==0]=0

plt.clf()
plt.plot(points[:,0],points[:,1],'*')
plt.savefig(wdfigs+'test',dpi=700)

records=MaskableList(records)
roadType=[rec.record[6] for rec in records]

maskP=np.array([rType=='primary' for rType in roadType])
maskS=np.array([rType=='secondary' for rType in roadType])
maskT=np.array([rType=='tertiary' for rType in roadType])

records1=records[maskP]
records2=records[maskS]
records3=records[maskT]

