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

##############################################
# Gridded Malnutrition
##############################################
print 'Malnutrition' 
try:
	mal=np.load(wdvars+'malnutrition_prevalence_2000-2015.npy')
	imageMask2=np.load(wdvars+'imageMask2.npy')
	
	imageMask3=np.zeros(shape=(nyears,10,imageMask2.shape[1],imageMask2.shape[2]),dtype=bool)
	for y in range(16):
		for i in range(10):
			imageMask3[y,i,:,:]=imageMask2[y] # imageMask2=(lat,lon,3)
	
	latm=latsubsaharan
	lonm=lonsubsaharan
	pixelsize=0.04166666666650246
	
	mal[y,:,:]=maltmp
	
	gridm = np.zeros(shape=(len(latm),len(lonm),2))
	for x in range(len(latm)):
		gridm[x,:,0]=latm[x]
	for x in range(len(lonm)):
		gridm[:,x,1]=lonm[x]
	
	gridMid=createMidpointGrid(gridm,pixelsize) # lon lat

except:
	print 'try command failed: retrieving malnutrition'

	mal=np.zeros(shape=(16,len(latsubsaharan),len(lonsubsaharan)))
	imageMask2=np.zeros(shape=(mal.shape),dtype=bool)
	imageMask3=np.zeros(shape=(imageMask2.shape[0],imageMask2.shape[1],imageMask2.shape[2],10),dtype=bool)
	nations=np.load(wdvars+'ID_grid')
	nationsMask=np.ma.getmask(nations)
	
	year=1999
	for y in range(16):
		year+=1
		ds=gdal.Open(wddata+'malnutrition/IHME_AFRICA_CGF_2000_2015_WASTING_MEAN_'+str(year)+'_PREVALENCE_Y2018M02D28.TIF')
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
		
		maltmp=ds.ReadAsArray()
		maltmp[maltmp>1]=0.0308212 # The average of this weird pixel's neighbors
		
		##### Scale to Africa #####
		maltmp=maltmp[latm<np.amax(latsubsaharan)+1e-2]
		latm=latm[latm<np.amax(latsubsaharan)+1e-2]
		maltmp=maltmp[latm>np.amin(latsubsaharan)-1e-2]
		latm=latm[latm>np.amin(latsubsaharan)-1e-2]
		
		maltmp=maltmp[:,lonm<np.amax(lonsubsaharan)+1e-2]
		lonm=lonm[lonm<np.amax(lonsubsaharan)+1e-2]
		maltmp=maltmp[:,lonm>np.amin(lonsubsaharan)-1e-2]
		lonm=lonm[lonm>np.amin(lonsubsaharan)-1e-2]
		########## Mask ###########
		latm=latsubsaharan
		lonm=lonsubsaharan
	
		mal[y,:,:]=maltmp
		
		latm=latm[::-1]
		gridm = np.zeros(shape=(len(latm),len(lonm),2))
		for x in range(len(latm)):
			gridm[x,:,0]=latm[x]
		for x in range(len(lonm)):
			gridm[:,x,1]=lonm[x]
		latm=latm[::-1]
		
		gridMid=createMidpointGrid(gridm,pixelsize) # lon lat
		
		imageMask2[y,mal[y]<0]=1
		imageMask2[y,africaMask1==1]=1
		imageMask2[y,900:1300,1450:]=1
		imageMask2[y,nationsMask==1]=1
		
		for i in range(10):
			imageMask3[y,i,:,:]=imageMask2[y] # imageMask2=(lat,lon,3)
		
		mal[y][mal[y]<0]=0

	np.save(wdvars+'malnutrition_prevalence_2000-2015.npy',mal)
	np.save(wdvars+'imageMask2.npy',imageMask2)
	
mal=np.ma.masked_array(mal,imageMask2)

if MakePlots:
	for y in range(nyears):
		plt.clf()
		plt.imshow(mal[y],cmap=cm.jet,vmax=0.3)
		plt.title(str(year)+': % Population with Severe Malnutrition')
		plt.xticks([])
		plt.yticks([])
		plt.colorbar()
		plt.savefig(wdfigs+'malnutrition_'+str(year),dpi=700)

	plt.clf()
	plt.imshow(np.mean(mal,axis=0),cmap=cm.jet,vmax=0.3)
	plt.title('Average % Population with Severe Malnutrition')
	plt.xticks([])
	plt.yticks([])
	plt.colorbar()
	plt.savefig(wdfigs+'malnutrition_avg',dpi=700)

##############################################
# Gridded Population
##############################################
print 'Population' 
try:
	pop5km1year=np.load(wdvars+'pop5km_2010.npy')

	pop5km=np.ma.zeros(shape=(nyears,len(latm),len(lonm)))
	for y in range(nyears):
		pop5km[y]=pop5km1year
except:
	print 'try command failed: retrieveing population'

	ds=gdal.Open(wddata+'population/gpw_v4_basic_demographic_characteristics_rev10_a000_004bt_2010_cntm_30_sec.tif')
	width1 = ds.RasterXSize
	height1 = ds.RasterYSize
	gt = ds.GetGeoTransform()
	minx = gt[0]
	miny = gt[3] + width1*gt[4] + height1*gt[5] 
	maxx = gt[0] + width1*gt[1] + height1*gt[2]
	maxy = gt[3] 
	pixelsizep=abs(gt[-1])
	
	# lat and lon
	latp=np.ones(shape=(height1))
	lonp=np.ones(shape=(width1))
	for w in range(width1):
		lonp[w]=minx+w*pixelsizep
	for h in range(height1):
		latp[h]=miny+h*pixelsizep
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
	###########################
	
	pop[pop<0]=0
	
	plt.clf()
	plt.imshow(pop,cmap=cm.gist_ncar_r,vmax=2000)
	plt.title('gridded population')
	plt.colorbar()
	plt.savefig(wdfigs+'pop5km',dpi=900)
	
	gridp=np.zeros(shape=(len(lonp),len(latp),2))
	for x in range(len(lonp)):
		gridp[x,:,0]=lonp[x]
	for x in range(len(latp)):
		gridp[:,x,1]=latp[y]
	
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
	np.save(wdvars+'pop5km_2010',pop5km)

pop5km=np.ma.masked_array(pop5km,imageMask2)

malnumber=mal*pop5km
malnumber=np.ma.masked_array(malnumber,imageMask2)

if MakePlots:
	plt.clf()
	plt.imshow(pop5km,cmap=cm.jet,vmax=2000)
	plt.title('gridded population')
	plt.colorbar()
	plt.savefig(wdfigs+'pop5km',dpi=700)

	for y in range(nyears):
		plt.clf()
		plt.imshow(malnumber,cmap=cm.nipy_spectral_r,vmax=500)
		plt.title('malnutrition number')
		plt.colorbar()
		plt.savefig(wdfigs+'malnumber',dpi=900)

######################################
# Countries
######################################
print 'Nations'
try:
	nations=np.load(wdvars+'ID_grid')
except:
	print 'try command failed: retrieveing nations'

	ds=gdal.Open(wddata+'boundaries/gpw-v4-national-identifier-grid-rev10_30_sec_tif/gpw_v4_national_identifier_grid_rev10_30_sec.tif')
	width = ds.RasterXSize
	height = ds.RasterYSize
	gt = ds.GetGeoTransform()
	minx = gt[0]
	miny = gt[3] + width*gt[4] + height*gt[5] 
	maxx = gt[0] + width*gt[1] + height*gt[2]
	maxy = gt[3] 
	pixelsizeNI=abs(gt[-1])
	
	latc=np.ones(shape=(height))
	lonc=np.ones(shape=(width))
	for w in range(width):
		lonc[w]=minx+w*pixelsizeNI
	for h in range(height):
		latc[h]=miny+h*pixelsizeNI
	
	latc=latc[::-1]
	
	nations=ds.ReadAsArray()
	##### Scale to Africa #####
	nations=nations[latc<np.amax(latm)+pixelsize]
	latc=latc[latc<np.amax(latm)+pixelsize]
	nations=nations[latc>np.amin(latm)]
	latc=latc[latc>np.amin(latm)]
	
	nations=nations[:,lonc<np.amax(lonm)+pixelsize]
	lonc=lonc[lonc<np.amax(lonm)+pixelsize]
	nations=nations[:,lonc>np.amin(lonm)]
	lonc=lonc[lonc>np.amin(lonm)]
	
	nations[nations==32767]=450
	##### Scale to 5km #####
	latl=np.radians(latc[::-1])+1.2
	lonl=np.radians(lonc)+1.2
	lut=RectSphereBivariateSpline(latl, lonl, nations)
	
	newLats,newLons=np.meshgrid(np.radians(latm[::-1])+1.2,np.radians(lonm)+1.2)
	nations=lut.ev(newLats.ravel(),newLons.ravel()).reshape((len(lonm),len(latm))).T
	
	nations=np.ma.masked_array(nations,imageMask2)
	nations=np.round(nations,0)
	nations.dump(wdvars+'ID_grid')

if MakePlots:
	plt.clf()
	plt.imshow(nations,cmap=cm.nipy_spectral,vmin=726,vmax=738)
	plt.yticks([])
	plt.xticks([])
	plt.title('Countries')
	plt.colorbar()
	plt.savefig(wdfigs +'countries',dpi=700)
	
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
	if len(np.where(nations==indexedcodes[i])[0])==0:
		print africanCountries[i]

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
try:
	malCountry = np.load(wdvars+'malCountry_2000-2015.npy')
except:
	print 'try command failed: calculating malCountry'
	malCountry=np.zeros(shape=(len(africanCountries),nyears))
	j=-1
	for i in indexedcodes:
		j+=1
		popWhole=np.sum(pop5km[nations==i])
		poptmp=np.ma.masked_array(pop5km,nations!=i)
		for y in range(nyears):
			maltmp=np.ma.masked_array(mal[y],nations!=i)
			malnumtmp=maltmp*poptmp
			malCountry[j,y]=np.sum(malnumtmp)/popWhole
	malCountry=malCountry*100
	
	malCountryGrid=np.zeros(shape=(nyears,len(latm),len(lonm)))
	j=-1
	for i in indexedcodes:
		j+=1
		for y in range(nyears):
			malCountryGrid[y,nations==i]=malCountry[j,y]
	
	malCountryGrid=np.ma.masked_array(malCountryGrid,imageMask2)
	np.save(wdvars+'malCountry_2000-2015',malCountry)

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

######################################
# Country Indices
######################################

################
index=0

f=open(wddata+'country_indices/gdp_per_cap.csv','r')

i=-2
g=-1
countries=[]
GDPperCap=np.zeros(shape=(len(africanCountries),nyears))
for line in f:
	i+=1
	if i==-1:
		continue
	line=line[:-2]
	line=line.replace('"','')
	tmp=np.array(line.split(','))
	country=tmp[0]
	if np.amax(country==np.array(africanCountries[:]))==0:
		continue
	countries.append(country)
	icountry=int(countryToi[country])
	g+=1
	j=43
	for y in range(nyears):
		j+=1
		try:
			GDPperCap[icountry,y]=float(tmp[j])
		except:
			exit()
Corr=np.zeros(shape=8)
Corr[index] = corr(np.ma.compressed(GDPperCap),np.ma.compressed(malCountry))

################
index+=1

f=open(wddata+'country_indices/percent_girls_educated.csv','r')

i=-2
g=-1
countries=[]
girlEd=np.zeros(shape=(len(africanCountries),nyears))
girlEdMask=np.zeros(shape=(len(africanCountries),nyears),dtype=bool)
for line in f:
	i+=1
	if i==-1:
		continue
	line=line[:-2]
	line=line.replace('"','')
	tmp=np.array(line.split(','))
	country=tmp[0]
	if np.amax(country==np.array(africanCountries[:]))==0:
		continue
	countries.append(country)
	icountry=int(countryToi[country])
	g+=1
	j=43
	for y in range(nyears):
		j+=1
		try:
			girlEd[icountry,y]=float(tmp[j])
		except:
			girlEdMask[icountry,y]=1

for i in range(len(girlEd[:,0])):
	edtmp=np.ma.masked_array(girlEd[i],girlEdMask[i])
	for y in range(nyears):
		ed1=0
		ed2=0
		if girlEd[i,y]==0:
			j1=y-1
			if girlEdMask[i,j1]==0:
				ed1=girlEd[i,j1]
			for j2 in range(y,nyears):
				if girlEdMask[i,j2]==0:
					ed2=girlEd[i,j2]
					break

			if ed1!=0 and ed2!=0:
				diff=ed2-ed1
				missingyears=j2-j1-1
				for k in range(y,j2):
					girlEd[i,k]=diff/missingyears+girlEd[i,k-1]
					girlEdMask[i,k]=0
			elif ed1!=0:
				girlEd[i,y]=girlEd[i,j1]
				girlEdMask[i,y]=0
			elif ed2!=0:
				girlEd[i,y]=girlEd[i,j2]
				girlEdMask[i,y]=0

girlEd=np.ma.masked_array(girlEd,girlEdMask)
malCountrytmp=np.ma.masked_array(malCountry,girlEdMask)
Corr[index] = corr(np.ma.compressed(girlEd),np.ma.compressed(malCountrytmp))

################
index+=1

f=open(wddata+'country_indices/fertility.csv','r')

i=-2
g=-1
countries=[]
fertility=np.zeros(shape=(len(africanCountries),nyears))
fertilityMask=np.zeros(shape=(len(africanCountries),nyears),dtype=bool)
for line in f:
	i+=1
	if i==-1:
		continue
	line=line[:-2]
	line=line.replace('"','')
	tmp=np.array(line.split(','))
	country=tmp[0]
	if np.amax(country==np.array(africanCountries[:]))==0:
		continue
	countries.append(country)
	icountry=int(countryToi[country])
	g+=1
	j=43
	for y in range(nyears):
		j+=1
		try:
			fertility[icountry,y]=float(tmp[j])
		except:
			fertilityMask[icountry,y]=1

Corr[index] = corr(np.ma.compressed(fertility),np.ma.compressed(malCountry))

################
index+=1

f=open(wddata+'country_indices/electricity.csv','r')

i=-2
g=-1
countries=[]
electricity=np.zeros(shape=(len(africanCountries),nyears))
electricityMask=np.zeros(shape=(len(africanCountries),nyears),dtype=bool)
for line in f:
	i+=1
	if i==-1:
		continue
	line=line[:-2]
	line=line.replace('"','')
	tmp=np.array(line.split(','))
	country=tmp[0]
	if np.amax(country==np.array(africanCountries[:]))==0:
		continue
	countries.append(country)
	icountry=int(countryToi[country])
	g+=1
	j=43
	for y in range(nyears):
		j+=1
		try:
			electricity[icountry,y]=float(tmp[j])
		except:
			electricityMask[icountry,y]=1

Corr[index] = corr(np.ma.compressed(electricity),np.ma.compressed(malCountry))

################
index+=1

f=open(wddata+'country_indices/cereal_yield.csv','r')

i=-2
g=-1
countries=[]
Yield=np.zeros(shape=(len(africanCountries),nyears))
yieldMask=np.zeros(shape=(len(africanCountries),nyears),dtype=bool)
for line in f:
	i+=1
	if i==-1:
		continue
	line=line[:-2]
	line=line.replace('"','')
	tmp=np.array(line.split(','))
	country=tmp[0]
	if np.amax(country==np.array(africanCountries[:]))==0:
		continue
	countries.append(country)
	icountry=int(countryToi[country])
	g+=1
	j=43
	for y in range(nyears):
		j+=1
		try:
			Yield[icountry,y]=float(tmp[j])
		except:
			yieldMask[icountry,y]=1

Yield=np.ma.masked_array(Yield,yieldMask)
malCountrytmp=np.ma.masked_array(malCountry,yieldMask)
Corr[index] = corr(np.ma.compressed(Yield),np.ma.compressed(malCountrytmp))

################
index+=1

f=open(wddata+'country_indices/under5mortality.csv','r')

i=-2
g=-1
countries=[]
mortality5=np.zeros(shape=(len(africanCountries),nyears))
mortality5Mask=np.zeros(shape=(len(africanCountries),nyears),dtype=bool)
for line in f:
	i+=1
	if i==-1:
		continue
	line=line[:-2]
	line=line.replace('"','')
	tmp=np.array(line.split(','))
	country=tmp[0]
	if np.amax(country==np.array(africanCountries[:]))==0:
		continue
	countries.append(country)
	icountry=int(countryToi[country])
	g+=1
	j=43
	for y in range(nyears):
		j+=1
		try:
			mortality5[icountry,y]=float(tmp[j])
		except:
			mortality5Mask[icountry,y]=1

Corr[index] = corr(np.ma.compressed(mortality5),np.ma.compressed(malCountry))

################
index+=1

f=open(wddata+'country_indices/iq.csv','r')

i=-2
g=-1
countries=[]
iq=np.zeros(shape=(len(africanCountries),nyears))
iqMask=np.zeros(shape=(len(africanCountries),nyears),dtype=bool)
for line in f:
	i+=1
	if i==-1:
		continue
	line=line[:-2]
	tmp=np.array(line.split(','))
	country=tmp[0]
	if np.amax(country==np.array(africanCountries[:]))==0:
		continue
	countries.append(country)
	icountry=int(countryToi[country])
	g+=1
	j=43
	for y in range(nyears):
		j+=1
		try:
			iq[icountry,y]=float(tmp[2])
		except:
			iqMask[icountry,y]=1

iq=np.ma.masked_array(iq,iqMask)
malCountrytmp=np.ma.masked_array(malCountry,iqMask)
Corr[index] = corr(np.ma.compressed(iq),np.ma.compressed(malCountry))

################
index+=1

f=open(wddata+'country_indices/corruption.csv','r')

i=-2
g=-1
countries=[]
corruption=np.zeros(shape=(len(africanCountries),nyears))
corruptionMask=np.zeros(shape=(len(africanCountries),nyears),dtype=bool)
for line in f:
    i+=1
    if i==-1:
        continue
    line=line[:-2]
    line=line.replace('"','')
    tmp=np.array(line.split(','))
    country=tmp[0]
    if np.amax(country==np.array(africanCountries[:]))==0:
        continue
    countries.append(country)
    icountry=int(countryToi[country])
    g+=1
    j=43
    for y in range(nyears):
        j+=1
        try:
            corruption[icountry,y]=float(tmp[6])
        except:
            corruptionMask[icountry,y]=1

Corr[index] = corr(np.ma.compressed(corruption),np.ma.compressed(malCountry))

################

#indices=[GDPperCap,girlEd,fertility,electricity,Yield,mortality5,iq,corruption]
indices=[girlEd,electricity,Yield,iq]
indexNames=['girlEd','electricity','Yield','iq']
xMulti=np.ma.zeros(shape=(len(africanCountries),nyears,len(indices)))

Mask = yieldMask + girlEdMask 
ydata=np.ma.compressed(np.ma.masked_array(malCountry,Mask))

xMultiC=np.zeros(shape=(len(ydata),len(indices)))
for i in range(len(indices)):
	xMultiC[:,i]=np.ma.compressed(np.ma.masked_array(indices[i],Mask))

malPred=np.zeros(shape=(10,141))
malActual=np.zeros(shape=(10,141))
for i in range(10):
	Xtrain,Xtest,Ytrain,Ytest=sklearn.model_selection.train_test_split(xMultiC,ydata,test_size=.2,train_size=.8)
	
	clf=linear_model.LinearRegression()
	clf.fit(Xtrain,Ytrain)
	clf=RandomForestRegressor()
	clf.fit(Xtrain,Ytrain)
	
	malPred[i]=clf.predict(Xtest)
	malActual[i]=Ytest

MultiCorr=corr(np.ma.compressed(malPred),np.ma.compressed(malActual))
print MultiCorr

plt.clf()
plt.plot(np.ma.compressed(malActual),np.ma.compressed(malPred),'b*')
plt.title('Acuracy of Malnutrition Predictions, Corr = '+str(round(MultiCorr,2)))
plt.xlabel('Actual Prevalence')
plt.ylabel('Predicted Prevalence')
plt.gca().set_aspect('equal', adjustable='box')
plt.grid(True)
plt.xlim([0,27])
plt.ylim([0,27])
plt.savefig(wdfigs+'malnutrition_predictions_accuracy.pdf')

for j in range(len(indices)):
	try:
		vars()[indexNames[j]+'Grid']=np.load(wdvars+indexNames[j]+'Grid_2000-2015')
	except:
		print 'try command failed: calculating',indexNames[j]+'Grid'

		vars()[indexNames[j]+'Grid']=np.zeros(shape=(nyears,len(latm),len(lonm)))
		k=-1
		for i in indexedcodes:
			k+=1
			for y in range(nyears):
				vars()[indexNames[j]+'Grid'][y,nations==i]=indices[j][k,y]
		
		vars()[indexNames[j]+'Grid'] = np.ma.masked_array(vars()[indexNames[j]+'Grid'], imageMask2)
			
		if MakePlots:
			plt.clf()
			plt.imshow(vars()[indexNames[j]+'Grid'][15],cmap=cm.jet_r)
			plt.colorbar()
			plt.yticks([])
			plt.xticks([])
			plt.title('2015 '+indexNames[j])
			plt.savefig(wdfigs+indexNames[j]+'.png',dpi=700)
		
		vars()[indexNames[j]+'Grid'].dump(wdvars+indexNames[j]+'Grid_2000-2015')

###########################################

###########################################
# Conflicts
###########################################
print 'Conflicts'
try:
	MPconflicts=np.load(wdvars+'MPconficts.npy')
except:
	conflictsAll=np.load(wdvars+'africa_fatalities')
	conflictCoordAll=np.load(wdvars+'conflictCoord')
	MPconflicts=np.zeros(shape=(nyears,len(latm),len(lonm)))
	
	year=1999
	k=-1
	for y in range(100,116):
		k+=1
		year+=1
		print '\n\n',year,'\n\n\n'
		conflicts=np.ma.compressed(conflictsAll[y:y+1])
		conflictCoord=np.zeros(shape=(len(conflicts),2))
		conflictCoord[:,0]=np.ma.compressed(conflictCoordAll[y:y+1,:,:,0])
		conflictCoord[:,1]=np.ma.compressed(conflictCoordAll[y:y+1,:,:,1])
	
		MPconflicts[k]=marketPotentials(conflicts,conflictCoord,gridMid,imageMask2[0,:,:])
	
	np.save(wdvars+'MPconficts',MPconflicts)
	
	year=1999
	k=-1
	for y in range(100,116):
		k+=1
		year+=1
		
		plt.clf()
		plt.imshow(MPconflicts[k],vmax=10000,cmap=cm.gist_heat_r)
		plt.colorbar()
		plt.xticks([])
		plt.yticks([])
		plt.title(str(year)+': Conflicts index')
		plt.savefig(wdfigs+'MPconflicts_'+str(year)+'.png',dpi=700)

MPconflicts=np.ma.masked_array(MPconflicts,MPconflicts==0)

######################################################
# Travel Time
######################################################
print 'Travel Time'
try:
	travel1year=np.load(wdvars+'travelAfrica')
	
	travel=np.ma.zeros(shape=(nyears,len(latm),len(lonm)))
	for y in range(nyears):
		travel[y]=travel1year
except: 
	#rRoads = shapefile.Reader(wddata+'openstreetmap/'+country+'/openstreetmap/nga_trs_roads_osm.shp')
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
	lnt=np.ones(shape=(width))
	for w in range(width):
		lont[w]=minx+w*pixelsizel
	for h in range(height):
		latt[h]=miny+h*pixelsizel
	
	latt=latt[::-1]
	if latm[0]<latm[-1]:
		latm=latm[::-1]
	
	travelWorld=ds.ReadAsArray()
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
	latl=np.radians(latt[::-1])+1.2
	lonl=np.radians(lont)+1.2
	lut=RectSphereBivariateSpline(latl, lonl, travel)
	
	newLats,newLons=np.meshgrid(np.radians(latm[::-1])+1.2,np.radians(lonm)+1.2)
	travel=lut.ev(newLats.ravel(),newLons.ravel()).reshape((len(lonm),len(latm))).T
	
	travel=np.ma.masked_array(travel,imageMask2)
	travel.dump(wdvars+'travelAfrica')

plt.clf()
plt.imshow(travel[0],cmap=cm.gist_ncar_r,vmax=400)
plt.title('Travel Time')
plt.colorbar()
plt.savefig(wdfigs+'travel',dpi=700)

###########################################
# Coast Lines
###########################################
print 'Coasts'
try:
	distToCoasts1year=np.load(wdvars+'Africa_distToCoasts.npy')
	coastLines=np.load(wdvars+'Africa_coastLines.npy')
	
	distToCoasts=np.ma.zeros(shape=(nyears,2,len(latm),len(lonm)))
	for y in range(nyears):
		for i in range(2):
			distToCoasts[y,i]=distToCoasts1year[:,:,i]
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
	
	np.save(wdvars+'Africa_coastLines.npy',coastLines)
	np.save(wdvars+'Africa_distToCoasts.npy',distToCoasts)

distToCoasts=np.ma.masked_array(distToCoasts,imageMask3[:,:2,:,:])
coastLines=np.ma.masked_array(coastLines,imageMask3[0,:2,:,:])

if MakePlots:
	plt.clf()
	plt.xticks([])
	plt.imshow(np.clip(distToCoasts[0,0,:,:],0,200))
	plt.colorbar()
	plt.xticks([])
	plt.yticks([])
	plt.title('Distance to Coastlines, km')
	plt.savefig(wdfigs+'distToCoasts'+str(1),dpi=700)


###########################################
# Land Classification
###########################################
'''
print 'Land Classification'
#try:
#	landCover=np.load(wdvars+'africa_landcover.npy')
#except:
ds=gdal.Open(wddata+'nigeria_landtype/MODIS_land_cover/allAfrica/2016_01_01.LC_Type3.tif')
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
latl=latl[::-1]

landCoversmall=ds.ReadAsArray()
landCoversmall=np.array(landCoversmall,dtype=float)
landCoversmall[landCoversmall==255]=0
#landCoversmall=np.ma.masked_array(landCoversmall,landCoversmall==255)
#landCoversmall=np.ma.masked_array(landCoversmall,landCoversmall==0)

#landCoversmall1=landCoversmall.filled(fill_value=np.NaN)
#landCoversmallFilled = replace_nans(landCoversmall1, 10, 0.5, 2, method='localmean')

##### Scale to Africa #####
landCoversmall=landCoversmall[latl<np.amax(latm)+1e-3]
latl=latl[latl<np.amax(latm)+1e-3]
landCoversmall=landCoversmall[latl>np.amin(latm)-1e-3]
latl=latl[latl>np.amin(latm)-1e-3]

landCoversmall=landCoversmall[:,lonl<np.amax(lonm)+pixelsize]
lonl=lonl[lonl<np.amax(lonm)+pixelsize]
landCoversmall=landCoversmall[:,lonl>np.amin(lonm)]
lonl=lonl[lonl>np.amin(lonm)]
	
gridl=np.zeros(shape=(len(latl),len(lonl),2))
for x in range(len(latl)):
	gridl[x,:,0]=latl[x]
for x in range(len(lonl)):
	gridl[:,x,1]=lonl[y]
	
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


#np.save(wdvars+'africa_landcover',landCover)

# urban
urbanS=np.zeros(shape=(landCoversmall.shape))
urbanS[landCoversmall==10]=1
# forest
forestS=np.zeros(shape=(landCoversmall.shape))
forestS[landCoversmall==5]=1
forestS[landCoversmall==6]=1
forestS[landCoversmall==7]=1
forestS[landCoversmall==8]=1
# shrublands
shrublandsS=np.zeros(shape=(landCoversmall.shape))
shrublandsS[landCoversmall==9]=1
shrublandsS[landCoversmall==1]=1
shrublandsS[landCoversmall==2]=1
shrublandsS[landCoversmall==4]=1
# croplands
croplandsS=np.zeros(shape=(landCoversmall.shape))
croplandsS[landCoversmall==3]=1
landCoversmall=np.array(landCoversmall)
exit()
#
########### Convert to 5km
##oldLat=np.radians(gridl[:,:,0])+1.2
##oldLon=np.radians(gridl[:,:,1])+1.2
##newLats,newLons=np.radians(latm[::-1])+1.2,np.radians(lonm)+1.2
##landCover=griddata(gridl+1.2,landCoversmall,(newLats,newLons),method='nearest')
#
## croplands
latt=np.radians(latl[::-1])+1.2
lont=np.radians(lonl)+1.2
lut=RegularGridInterpolator((latt,lont), landCoversmall, 'nearest')

newLats,newLons=np.meshgrid(np.radians(latm[::-1])+1.2,np.radians(lonm)+1.2)
landCover=lut((newLats.ravel(),newLons.ravel())).reshape((len(lonm),len(latm))).T
croplands=np.round(landCover,0)
croplands=np.ma.masked_array(landCover,imageMask2)
###########
#
#
#
#plt.clf()
#plt.imshow(croplands,cmap=cm.binary)
#plt.colorbar()
##plt.xticks([])
##plt.yticks([])
#plt.title('Nigeria Land Cover')
#plt.savefig(wdfigs+'landCover',dpi=700)
exit()
'''

###########################################
# Machine learning 
###########################################
print 'Machine Learning'
#indices=[pop5km,distToCoasts[:,:,0],distToCoasts[:,:,1],MPconflicts,travel,GDPgrid,urban,forest,shrublands,croplands]
indices=[pop5km,distToCoasts[:,0,:,:],distToCoasts[:,1,:,:],MPconflicts,travel,girlEdGrid,electricityGrid,YieldGrid,iqGrid]
#indices=[pop5km,np.clip(MPconflicts,0,40000),distToCoasts[:,:,1]]
#indices=[distToCoasts[:,:,0]]
xMulti=np.ma.zeros(shape=(nyears,len(latm),len(lonm),len(indices)))
ydata=np.ma.compressed(mal)

for i in range(len(indices)):
	xMulti[:,:,:,i]=indices[i]

###########################
# Split into boxes 
###########################
if latm[0]>latm[-1]:
	latm=latm[::-1]
ysplit=len(latm)/20.
xsplit=len(lonm)/20.

yboxes=np.arange(0,len(latm)+1,ysplit)
xboxes=np.arange(0,len(lonm)+1,xsplit)

boxesGrid=np.zeros(shape=(gridm.shape))
for i in range(len(yboxes)-1):
	boxesGrid[:,:,0][gridm[:,:,0]>latm[int(np.round(yboxes[i],0))]]=i
for i in range(len(xboxes)-1):
	boxesGrid[:,:,1][gridm[:,:,1]>lonm[int(np.round(xboxes[i],0))]]=i

plt.clf()
plt.imshow(boxesGrid[:,:,0])
plt.title('Latitude Boxes')
plt.colorbar()
plt.savefig(wdfigs+'latboxes.pdf')
plt.clf()
plt.imshow(boxesGrid[:,:,1])
plt.title('Longitude Boxes')
plt.colorbar()
plt.savefig(wdfigs+'lonboxes.pdf')

boxes=np.zeros(shape=(boxesGrid[:,:,0].shape))
k=-1
for i in range(len(yboxes)):
	for j in range(len(xboxes)):
		k+=1
		ytrue=boxesGrid[:,:,0]==i
		xtrue=boxesGrid[:,:,1]==j
		true=ytrue*xtrue
		boxes[true]=k

boxes=np.array(boxes)
boxNums=np.unique(boxes)
boxNums=np.array(boxNums,dtype=int)

boxMask=np.zeros(shape=(boxes.shape))
for i in boxNums:
	masktmp=np.zeros(shape=(boxes.shape),dtype=bool)
	masktmp[boxes==i]=1
	if np.sum(1-imageMask2[masktmp])<100:
		boxMask[masktmp==1]=1
boxes=np.ma.masked_array(boxes,boxMask)

boxNums=np.unique(boxes)

trainNums,testNums=sklearn.model_selection.train_test_split(boxNums,test_size=.2,train_size=.8)
trainPixelsAll=[boxes==np.array(trainNums[i]) for i in range(len(trainNums))]
trainPixels=np.sum(trainPixelsAll,axis=0)
testPixelsAll=[boxes==np.array(testNums)[i] for i in range(len(testNums))] 
testPixels=np.sum(testPixelsAll,axis=0) 
trainPixels=1-trainPixels
testPixels=1-testPixels
exit()

maltrain=np.ma.array(mal)
maltrain=np.ma.masked_array(maltrain,trainPixels)
maltest=np.ma.array(mal)
maltest=np.ma.masked_array(maltest,testPixels)

plt.clf()
plt.imshow(maltrain,cmap=cm.jet,vmin=0,vmax=0.3)
plt.title('Test Pixels')
plt.colorbar()
plt.savefig(wdfigs+'maltrain.png',dpi=700)
plt.clf()
plt.imshow(maltest,cmap=cm.jet,vmin=0,vmax=0.3)
plt.title('Test Pixels')
plt.colorbar()
plt.savefig(wdfigs+'maltest.png',dpi=700)

boxMask=np.ma.masked_array(boxMask,testPixels)
plt.clf()
plt.imshow(boxMask,cmap=cm.seismic)
plt.title('Boxes')
plt.colorbar()
plt.savefig(wdfigs+'trainboxes.png',dpi=700)

plt.clf()
plt.imshow(boxes,cmap=cm.seismic)
plt.title('Boxes')
plt.colorbar()
plt.savefig(wdfigs+'boxes.png',dpi=700)
exit()

trainMask=np.ma.getmask(maltrain)
testMask=np.ma.getmask(maltest)
trainMask3=np.zeros(shape=(len(latm),len(lonm),len(indices)),dtype=bool)
testMask3=np.zeros(shape=(len(latm),len(lonm),len(indices)),dtype=bool)
for i in range(len(xMulti[0,0,:])):
	trainMask3[:,:,i]=trainMask
	testMask3[:,:,i]=testMask

xMultiTest=np.ma.masked_array(xMulti, trainMask3)
xMultiTrain=np.ma.masked_array(xMulti, testMask3)

xMultiTrainC=np.zeros(shape=(len(np.ma.compressed(xMultiTrain[:,:,0])),len(xMulti[0,0,:])))
xMultiTestC=np.zeros(shape=(len(np.ma.compressed(xMultiTest[:,:,0])),len(xMulti[0,0,:])))
for i in range(len(xMulti[0,0,:])):
	xMultiTrainC[:,i]=np.ma.compressed(xMultiTrain[:,:,i])
	xMultiTestC[:,i]=np.ma.compressed(xMultiTest[:,:,i])
#ydataTrain=np.ma.compressed(maltrain)
ydataTrain=np.ma.compressed(maltest)
#ydataTest=np.ma.compressed(maltest)
ydataTest=np.ma.compressed(maltrain)

clf=RandomForestRegressor()
clf.fit(xMultiTrainC,ydataTrain)
malPred=clf.predict(xMultiTestC)

#clf=linear_model.LinearRegression()
#clf.fit(xMultiTrain,ydataTrain)
#
#malPred=clf.predict(xMultiTest)
Corr=corr(malPred,ydataTest)
print Corr

x,y=ydataTest*100,malPred*100
slope,b=np.polyfit(x,y,1)
yfit=slope*x+b
error=np.mean((y-x)/x)*100
heatmap, xedges, yedges = np.histogram2d(x, y, bins=(200,200),normed=True)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
heatmap = heatmap.T

plt.clf()
fig = plt.figure(figsize=(7, 5))
ax = fig.add_subplot(111, title='Accuracy of Predicted Malnutrition, Corr = '+str(round(Corr,2))+', Avg Error ='+str(round(error,2))+'%')
X, Y = np.meshgrid(xedges, yedges)
ax.pcolormesh(X, Y, heatmap,cmap=cm.terrain_r)
plt.xlabel('Actual Malnutrition, % population')
plt.ylabel('Predicted Malnutrition, % population')
ax.set_xlim([0,30])
ax.set_ylim([0,30])
plt.savefig(wdfigs+'accuracy_predictions_heatmap.pdf') 

exit()

malPred=np.zeros(shape=(MPconflicts.shape))
malPred=np.ma.masked_array(malPred,imageMask2)
testtmp=imageMask2==testMask
traintmp=imageMask2==trainMask
x1test=0
x1train=0
for i in range(len(malPred[:,0])):
	x2test=np.sum(testtmp[i]==False)+x1test
	x2train=np.sum(traintmp[i]==False)+x1train
	if x2train==x1train:
		continue
	malPred[i][traintmp[i]==False]=ydataTrain[x1train:x2train]
	if x1test!=x2test:
		malPred[i][testtmp[i]==False]=clf.predict(xMultiTestC[x1test:x2test])
	x1train=x2train
	x1test=x2test
	

plt.clf()
plt.imshow(malPred,cmap=cm.jet,vmin=0,vmax=0.3)
plt.colorbar()
plt.yticks([])
plt.xticks([])
plt.title('Predicted Malnutrition, Percent of Population')
plt.savefig(wdfigs+'malPred_Africa',dpi=700)
Corr=corr(np.ma.compressed(malPred),np.ma.compressed(mal))


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
travel=ds.ReadAsArray()
travel[travel<0]=-10

africaMask1=np.zeros(shape=(travel.shape),dtype=int)
africaMask1[travel<0]=1

######### 5km #########
lattmp=latm[latm<np.amax(latc)+pixelsize]
lattmp=lattmp[lattmp>np.amin(latc)]
travel=travel[latc>latm[-1]]
africaMask1=africaMask1[latc>latm[-1]]
latc=latc[latc>latm[-1]]

latl=np.radians(latc[::-1])+1.2
lonl=np.radians(lonc)+1.2
lut=RectSphereBivariateSpline(latl, lonl, travel)
newLats,newLons=np.meshgrid(np.radians(lattmp[::-1])+1.2,np.radians(lonm)+1.2)
travel=lut.ev(newLats.ravel(),newLons.ravel()).reshape((len(lonm),len(lattmp))).T

lut=RectSphereBivariateSpline(latl, lonl, africaMask1)
newLats,newLons=np.meshgrid(np.radians(lattmp[::-1])+1.2,np.radians(lonm)+1.2)
africaMask1=lut.ev(newLats.ravel(),newLons.ravel()).reshape((len(lonm),len(lattmp))).T
#######################
africaMask1=np.round(africaMask1,0)

travel=np.ma.masked_array(travel,africaMask1)

plt.clf()
plt.imshow(travel,cmap=cm.gist_ncar_r,vmin=0,vmax=30)
plt.title('Travel Time')
plt.yticks([])
plt.xticks([])
plt.colorbar()
plt.savefig(wdfigs+'travel',dpi=700)
