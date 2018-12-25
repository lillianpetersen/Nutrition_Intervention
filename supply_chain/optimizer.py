from sys import exit
import csv
#import supplychainpy
from math import sqrt
from sys import exit
import pickle
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import matplotlib.mlab as mlab
#from mpl_toolkits.basemap import Basemap
import os
from scipy.stats import norm
import matplotlib as mpl
from matplotlib.patches import Polygon
import matplotlib.cm as cm
import random
from libtiff import TIFF
from osgeo import gdal
from scipy.interpolate import RectSphereBivariateSpline

def createMidpointGrid(grid,pixelsize):
	gridMid=np.zeros(shape=(grid.shape))
	for ilon in range(len(grid[:,0,0])):
		gridMid[ilon,:,0]=np.mean([grid[ilon,0,0],grid[ilon,0,0]+pixelsize])
	for ilat in range(len(grid[0,:,1])):
		gridMid[:,ilat,1]=np.mean([grid[:,ilat,1],grid[:,ilat,1]+pixelsize])
	return gridMid


wddata='C:/Users/garyk/Documents/code/riskAssessmentFromPovertyEstimations/supply_chain/data/'
wdvars='C:/Users/garyk/Documents/code/riskAssessmentFromPovertyEstimations/supply_chain/variables/'
wdfigs='C:/Users/garyk/Documents/code/riskAssessmentFromPovertyEstimations/supply_chain/figures/'

##population
wddata='C:/Users/garyk/Documents/code/riskAssessmentFromPovertyEstimations/data/'
# 
pop20Tif=TIFF.open(wddata+'GHS_SMOD_POP2015HDC_GLOBE_R2016A_54009_1k_v1_0\GHS_SMOD_POP2015HDC_GLOBE_R2016A_54009_1k_v1_0/GHS_SMOD_POP2015HDC_GLOBE_R2016A_54009_1k_v1_0.tif',mode='r')
mc=gdal.Open(wddata+'GHS_SMOD_POP2015HDC_GLOBE_R2016A_54009_1k_v1_0\GHS_SMOD_POP2015HDC_GLOBE_R2016A_54009_1k_v1_0/GHS_SMOD_POP2015HDC_GLOBE_R2016A_54009_1k_v1_0.tif')
# pop20Tif=TIFF.open('C:/Users/garyk/Documents/python_code/supply_chain/data/AFRIPOP_2015_U5/afripop2015u5.tif',mode='r')
# mc=gdal.Open('C:/Users/garyk/Documents/python_code/supply_chain/data/AFRIPOP_2015_U5/afripop2015u5.tif')
# pop20Tif = TIFF.open(wddata + 'AFR_PPP_A0004_M_2020_adj_v5.tif', mode='r')
# mc = gdal.Open(wddata + 'AFR_PPP_A0004_M_2020_adj_v5.tif')
width = mc.RasterXSize
height = mc.RasterYSize
plt.clf()
pop20=pop20Tif.read_image()
imageMaskPo = pop20<0
pop20=np.ma.masked_array(pop20,imageMaskPo)
plt.imshow(pop20,cmap=cm.binary,vmin=0,vmax=1)
plt.yticks([])
plt.xticks([])
plt.title('Children under 4 in 2020')
plt.colorbar()
plt.savefig(wdfigs +'2020MChildPop',dpi=700)

exit()
################################maybe for the bivariate should i translate it and translate it back?

###population
# pop20Tif = TIFF.open(wddata + 'IHME_AFRICA_CGF_2000_2015_WASTING_MEAN_2015_PREVALENCE_Y2018M02D28.tif', mode='r')
# mc = gdal.Open(wddata + 'IHME_AFRICA_CGF_2000_2015_WASTING_MEAN_2015_PREVALENCE_Y2018M02D28.tif')
# width = mc.RasterXSize
# height = mc.RasterYSize
# plt.clf()
# pop20=pop20Tif.read_image()*100
# imageMaskPo = pop20<0
# pop20=np.ma.masked_array(pop20,imageMaskPo)
# plt.imshow(pop20,cmap=cm.jet_r,vmin=0,vmax=100)
# plt.yticks([])
# plt.xticks([])
# plt.title('Wasting in 2015')
# plt.colorbar()
# plt.savefig(wdfigs +'Wasting_2015',dpi=700)

###stunting
stuntingTif = TIFF.open(wddata + 'IHME_AFRICA_CGF_2000_2015_STUNTING_MEAN_2015_PREVALENCE_Y2018M02D28.tif', mode='r')
mc = gdal.Open(wddata + 'IHME_AFRICA_CGF_2000_2015_STUNTING_MEAN_2015_PREVALENCE_Y2018M02D28.tif')
width1 = mc.RasterXSize
height1 = mc.RasterYSize
gt = mc.GetGeoTransform()
minx = gt[0]
miny = gt[3] + width1*gt[4] + height1*gt[5] 
maxx = gt[0] + width1*gt[1] + height1*gt[2]
maxy = gt[3] 
pixelsize=abs(gt[-1])

lat=np.ones(shape=(height1))
lon=np.ones(shape=(width1))
for w in range(width1):
	lon[w]=minx+w*pixelsize
for h in range(height1):
	lat[h]=miny+h*pixelsize
lat=lat[::-1]

stunting15 = stuntingTif.read_image()*100
imageMask = stunting15<0

stunting15=np.ma.masked_array(stunting15,imageMask)
stunting15=stunting15[~imageMask.all(axis=1)]
stunting15=stunting15[:,~imageMask.all(axis=0)]

# plt.clf()
# plt.imshow(stunting15,cmap=cm.jet_r,vmin=0,vmax=50)
# plt.yticks([])
# plt.xticks([])
# plt.title('Stunting Prevalence in 2015')
# plt.colorbar()
# plt.savefig(wdfigs +'stunting_2015',dpi=700)
###pops

mchildTif = TIFF.open(wddata + 'AFR_PPP_A0004_M_2015_adj_v5.tif')
mchild=mchildTif.read_image()
mc = gdal.Open(wddata + 'AFR_PPP_A0004_M_2015_adj_v5.tif')
width = mc.RasterXSize
height = mc.RasterYSize
gt = mc.GetGeoTransform()
minx = gt[0]
miny = gt[3] + width*gt[4] + height*gt[5] 
maxx = gt[0] + width*gt[1] + height*gt[2]
maxy = gt[3] 
pixelsize=abs(gt[-1])

latP=np.ones(shape=(height))
lonP=np.ones(shape=(width))
for w in range(width):
	lonP[w]=minx+w*pixelsize
for h in range(height):
	latP[h]=miny+h*pixelsize
latP=latP[::-1]

mchildTif = TIFF.open(wddata + 'AFR_PPP_A0004_M_2015_adj_v5.tif')
mchild=mchildTif.read_image()
imageMask=mchild<0

mchild=np.ma.masked_array(mchild,imageMask)
mchild=mchild[~imageMask.all(axis=1)]
mchild=mchild[:,~imageMask.all(axis=0)]

lonM=np.ma.masked_array(lonP,imageMask.all(axis=0))
lonM=np.ma.compressed(lonM)
latM=np.ma.masked_array(latP,imageMask.all(axis=1))
latM=np.ma.compressed(latM)

## Masks
imageMask2=imageMask # imageMask2=(lat,lon)
imageMask2=imageMask2[~imageMask2.all(axis=1)]
imageMask2=imageMask2[:,~imageMask2.all(axis=0)]
imageMask1=np.swapaxes(imageMask2,0,1) # imageMask1=(lon,lat)
imageMask3=np.zeros(shape=(imageMask2.shape[0],imageMask2.shape[1],10),dtype=bool)
for i in range(10):
	imageMask3[:,:,i]=imageMask2 # imageMask2=(lat,lon,3
	
## grid
grid = np.zeros(shape=(len(lonM),len(latM),2))
for x in range(len(lonM)):
	grid[x,:,0]=lonM[x]
for y in range(len(latM)):
	grid[:,y,1]=latM[y]

########################
# Final Downscaled Grid
########################
newLats,newLons=lat,lon
pixelsizeCoarse=pixelsize
latMr=latM[::-1]

gridCoarse=np.zeros(shape=(len(newLons),len(newLats),2))
for x in range(len(newLons)):
	gridCoarse[x,:,0]=newLons[x]
for y in range(len(newLats)):
	gridCoarse[:,y,1]=newLats[y]

## Plot old and new grid
# plt.clf()
# m=5*3
# plt.plot(np.ma.compressed(grid[:m,:m,0]),np.ma.compressed(grid[:m,:m,1]),'*k')
# plt.plot(np.ma.compressed(gridCoarse[:3,:3,0]),np.ma.compressed(gridCoarse[:3,:3,1]),'*r')
# plt.savefig(wdfigs + 'old_new_grid.pdf')
##

## Find New Mask
oldLat,oldLon=np.radians(latM[::-1]),np.radians(lonM)
oldLat2=oldLat[oldLat>0]
oldLon2=oldLon[oldLon>0]
imageMask2=imageMask2[oldLat>0]
imageMask2=imageMask2[:,oldLon>0]
mchild1=mchild[oldLat>0]
mchild1=mchild1[:,oldLon>0]


#oldLon=oldLon
#oldLat=oldLat

lutMask=RectSphereBivariateSpline(oldLat2, oldLon2, imageMask2)

lutChild=(oldLat2, oldLon2, mchild1)

newLats1,newLons1=np.radians(newLats[::-1]),np.radians(newLons)
#newLons=newLons+1.2
#newLats=newLats+1.2
newLats1=newLats1[newLats1>0]
newLons1=newLons1[newLons1>0]
newLats1m,newLons1m=np.meshgrid(newLats1,newLons1)

maskCoarse=lutMask.ev(newLats1m.ravel(),newLons1m.ravel()).reshape((len(newLons1),len(newLats1))).T
#maskCoarse=np.array(np.round(maskCoarse,0),dtype=bool)

mchildCoarse=lutChild.ev(newLats1m.ravel(),newLons1m.ravel()).reshape((len(newLons1),len(newLats1))).T
plt.clf()
plt.imshow(mchildCoarse,cmap=cm.coolwarm,vmin=0,vmax=1)
plt.yticks([])
plt.xticks([])
plt.title('Number of MChildren Age 0-4 2015')
plt.colorbar()
plt.savefig(wdfigs + 'Mchildren_test2',dpi=700)
plt.show()
exit()

mchildCoarse=np.ma.masked_array(mchildCoarse,maskCoarse)

plt.clf()
plt.imshow(mchildCoarse,cmap=cm.coolwarm,vmin=0,vmax=1)
plt.yticks([])
plt.xticks([])
plt.title('Number of MChildren Age 0-4 2015')
plt.colorbar()
plt.savefig(wdfigs + 'Mchildren_population',dpi=700)

exit()

fchildTif = TIFF.open(wddata + 'AFR_PPP_A0004_F_2015_adj_v5.tif')
fchild=mchildTif.read_image()

imageMask = mchild<0

mchild=np.ma.masked_array(mchild,imageMask)
mchild=mchild[~imageMask.all(axis=1)]
mchild=mchild[:,~imageMask.all(axis=0)]

fchild=np.ma.masked_array(fchild,imageMask)
fchild=fchild[~imageMask.all(axis=1)]
fchild=fchild[:,~imageMask.all(axis=0)]

child=mchild+fchild

plt.clf()
plt.imshow(child,cmap=cm.jet,vmin=0,vmax=40)
plt.yticks([])
plt.xticks([])
plt.title('Number of Children Age 0-4 2015')
plt.colorbar()
plt.savefig(wdfigs + 'children_population',dpi=700)

stuntchild = child*stunting15

plt.clf()
plt.imshow(child,cmap=cm.jet,vmin=0,vmax=30)
plt.yticks([])
plt.xticks([])
plt.title('Number of Stunted Children Age 0-4 2015')
plt.colorbar()
plt.savefig(wdfigs + 'stunted_children_population',dpi=700)